import numpy as np
from funcoes import MaxLoad


def pso(n_particulas, n_iteracoes, n_pgs):
    posicoes = np.random.uniform(Pgmin, Pgmax, (n_particulas, n_pgs))
    velocidades = np.zeros((n_particulas,n_pgs))
    fitness = np.array([MaxLoad(tipo, posicao) for posicao in posicoes])
    pbest = np.copy(posicoes)
    valor_pbest = np.copy(fitness)
    posicao_gbest = np.zeros(n_pgs)
    valor_gbest = -np.inf
    for i in range(n_particulas):
        if fitness[i] > valor_gbest:
            valor_gbest = fitness[i]
            posicao_gbest = np.copy(posicoes[i])
    print(f'Iteração 0 concluída! Carregamento máximo: {valor_gbest:.4f}')
    for i in range(n_iteracoes):
        r1 = np.random.uniform(0, 1, (n_particulas, n_pgs))
        r2 = np.random.uniform(0, 1, (n_particulas, n_pgs))
        velocidades = w * velocidades + c1 * r1 * (pbest - posicoes) + c2 * r2 * (posicao_gbest - posicoes)
        posicoes += velocidades
        for j in range(n_particulas):
            posicoes[j, :] = np.clip(posicoes[j, :], Pgmin, 9999)
        fitness = np.array([MaxLoad(tipo, posicao) for posicao in posicoes])
        for j in range(n_particulas):
            if fitness[j] > valor_pbest[j]:
                valor_pbest[j] = fitness[j]
                pbest[j] = np.copy(posicoes[j])
        melhor_fitness_idx = np.argmax(fitness)
        if fitness[melhor_fitness_idx] > valor_gbest:
            valor_gbest = fitness[melhor_fitness_idx]
            posicao_gbest = np.copy(posicoes[melhor_fitness_idx])
        print(f'Iteração {i + 1} concluída! Carregamento máximo: {valor_gbest:.4f}')
    return valor_gbest, posicao_gbest

tipo = 'ATIVO'
n_pgs = 1
Pgmin = 0
Pgmax = 10

w = 0.8
c1 = 1.25
c2 = 1.25

n_particulas = 10
n_iteracoes = 20

variaveis = ['Pg2']

carregamento_max, valor_vpg = pso(n_particulas, n_iteracoes, n_pgs)

print(valor_vpg)

print('Resultado final')
for i in range(len(variaveis)):
    print(f'{variaveis[i]} = {valor_vpg[i]:.4f} pu')
print(f'Carregamento máximo: {carregamento_max}')



