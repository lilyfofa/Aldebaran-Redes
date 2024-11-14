import numpy as np
from funcoes import MaxLoad, YBus
from time import time


def desvio(posicao):
    tensoes = posicao[:n_tensoes]
    referencia = 1.0
    desvio_medio = np.mean(np.abs(np.array(tensoes) - referencia))
    return desvio_medio


def pso(n_particulas, n_iteracoes, n_tensoes):
    ybus = YBus()
    posicoes = np.random.uniform(Vmin, Vmax, (n_particulas, n_tensoes))
    velocidades = np.zeros((n_particulas,n_tensoes))
    fitness = np.array([MaxLoad(tipo, posicao, ybus) for posicao in posicoes])
    pbest = np.copy(posicoes)
    valor_pbest = np.copy(fitness)
    posicao_gbest = np.zeros(n_tensoes)
    valor_gbest = -np.inf
    for i in range(n_particulas):
        if fitness[i] > valor_gbest:
            valor_gbest = fitness[i]
            posicao_gbest = np.copy(posicoes[i])
        elif fitness[i] == valor_gbest and desvio(posicoes[i]) < desvio(posicao_gbest):
            posicao_gbest = np.copy(posicoes[i])
    print(f'Iteração 0 concluída! Carregamento máximo: {valor_gbest:.4f}. Posição: {posicao_gbest}')
    for i in range(n_iteracoes):
        r1 = np.random.uniform(0, 1, (n_particulas, n_tensoes))
        r2 = np.random.uniform(0, 1, (n_particulas, n_tensoes))
        velocidades = w * velocidades + c1 * r1 * (pbest - posicoes) + c2 * r2 * (posicao_gbest - posicoes)
        posicoes += velocidades
        for j in range(n_particulas):
            posicoes[j, :] = np.clip(posicoes[j, :], Vmin, Vmax)
        fitness = np.array([MaxLoad(tipo, posicao, ybus) for posicao in posicoes])
        for j in range(n_particulas):
            if fitness[j] > valor_pbest[j]:
                valor_pbest[j] = fitness[j]
                pbest[j] = np.copy(posicoes[j])
        melhor_fitness_idx = np.argmax(fitness)
        if fitness[melhor_fitness_idx] > valor_gbest:
            valor_gbest = fitness[melhor_fitness_idx]
            posicao_gbest = np.copy(posicoes[melhor_fitness_idx])
        elif fitness[melhor_fitness_idx] == valor_gbest and desvio(posicoes[melhor_fitness_idx]) < desvio(posicao_gbest):
            posicao_gbest = np.copy(posicoes[melhor_fitness_idx])
        print(f'Iteração {i + 1} concluída! Carregamento máximo: {valor_gbest:.4f}. Posição: {posicao_gbest}')
    return valor_gbest, posicao_gbest

tipo = 'REATIVO'
n_tensoes = 5

Vmin = 0.95
Vmax = 1.10

w = 0.8
c1 = 1.25
c2 = 1.25

n_particulas = 10
n_iteracoes = 20
n_repeticoes = 1

variaveis = ['V1', 'V2', 'V3', 'V6', 'V8']

for i in range(n_repeticoes):
    t1 = time()
    carregamento_max, valor_vpg = pso(n_particulas, n_iteracoes, n_tensoes)
    t2 = time()

    print(f'Resultado da repetição {i+1}')
    for i in range(len(variaveis)):
        print(f'{variaveis[i]} = {valor_vpg[i]:.4f} pu')
    print(f'Carregamento máximo: {carregamento_max:.4f}')
    print(f'Desvio médio: {desvio(valor_vpg):.4f}')
    print(f'Tempo de execução: {t2-t1:.2f} s')
