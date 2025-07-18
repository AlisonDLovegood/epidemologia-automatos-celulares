import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Parâmetros CA-SEIR
sigma = 0.5    # incubação em 2 dias (1/2)
gamma = 0.14   # recuperação em 7 dias (1/7)
beta  = 0.18   # transmissão por vizinho (mediana_Influenza_sazonal × sigma == 1.28 × 0.14)

# Inicialização do grid
N = 100
grid = np.zeros((N, N), dtype=int)  # 0=S, 1=E, 2=I, 3=R
grid[N//2, N//2] = 2

def count_infected_neighbors(grid, x, y):
    N = grid.shape[0]
    count = 0
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            nx, ny = x + dx, y + dy
            if 0 <= nx < N and 0 <= ny < N and grid[nx, ny] == 2:
                count += 1
    return count

def update_seir_ca(grid, beta, sigma, gamma):
    N = grid.shape[0]
    new_grid = grid.copy()
    for x in range(N):
        for y in range(N):
            state = grid[x, y]
            if state == 0:  # S
                nI = count_infected_neighbors(grid, x, y)
                p_inf = 1 - (1 - beta)**nI
                if np.random.random() < p_inf:
                    new_grid[x, y] = 1
            elif state == 1:  # E
                if np.random.random() < sigma:
                    new_grid[x, y] = 2
            elif state == 2:  # I
                if np.random.random() < gamma:
                    new_grid[x, y] = 3
    return new_grid

def animate_ca(grid, beta, sigma, gamma, steps, pause_time):
    plt.figure(figsize=(7,7))
    for step in range(steps):
        plt.clf()
        plt.title(f"Passo {step}")
        plt.imshow(grid, cmap=ListedColormap(['white','gold','crimson','limegreen']), vmin=0, vmax=3)
        plt.pause(pause_time)
        grid = update_seir_ca(grid, beta, sigma, gamma)
    plt.show()

# Nova função para coletar dados e gerar gráfico de contagens
def simulate_seir_ca(steps):
    g = grid.copy()
    counts = []
    for _ in range(steps):
        unique, cnt = np.unique(g, return_counts=True)
        count_dict = {s: 0 for s in (0,1,2,3)}
        count_dict.update(dict(zip(unique, cnt)))
        counts.append((count_dict[0], count_dict[1], count_dict[2], count_dict[3]))
        g = update_seir_ca(g, beta, sigma, gamma)
    return np.array(counts), g

def plot_seir_counts(counts):
    plt.figure(figsize=(8,4))
    plt.plot(counts[:,0], label='Suscetível', color='blue')
    plt.plot(counts[:,1], label='Exposto',       color='orange')
    plt.plot(counts[:,2], label='Infectado',     color='red')
    plt.plot(counts[:,3], label='Recuperado',    color='green')
    plt.xlabel('Tempo (steps)')
    plt.ylabel('Número de células')
    plt.title('Evolução SEIR no CA')
    plt.legend()
    plt.grid(True)
    plt.show()

animate_ca(grid, beta, sigma, gamma, steps=250, pause_time=0.05)

counts, final_grid = simulate_seir_ca(250)
plot_seir_counts(counts)
