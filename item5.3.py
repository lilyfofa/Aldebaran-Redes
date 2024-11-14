import numpy as np
from funcoes import MaxLoad, YBus
from time import time


def performance(posicao):
    tensoes = posicao[:n_tensoes]
    potencia = posicao[n_tensoes:]
    ref = 1.0
    a = 0.8
    b = 0.2
    valor = a * np.mean(np.abs(np.array(tensoes)-ref)) + b * np.mean(potencia)
    return valor

def pso(n_particulas, n_iteracoes, n_pgs, n_tensoes):
    ybus = YBus()
    tensoes = np.random.uniform(Vmin, Vmax, (n_particulas, n_tensoes))
    geracoes = np.random.uniform(Pgmin, Pgmax, (n_particulas, n_pgs))
    posicoes = np.block([tensoes, geracoes])
    velocidades = np.zeros((n_particulas,n_pgs+n_tensoes))
    fitness = np.array([MaxLoad(tipo, posicao, ybus) for posicao in posicoes])
    pbest = np.copy(posicoes)
    valor_pbest = np.copy(fitness)
    posicao_gbest = np.zeros(n_pgs + n_tensoes)
    valor_gbest = -np.inf
    for i in range(n_particulas):
        if fitness[i] > valor_gbest:
            valor_gbest = fitness[i]
            posicao_gbest = np.copy(posicoes[i])
        elif fitness[i] == valor_gbest and performance(posicoes[i]) < performance(posicao_gbest):
            posicao_gbest = np.copy(posicoes[i])
    print(f'Iteração 0 concluída! Carregamento máximo: {valor_gbest:.4f}. Posição: {posicao_gbest}')
    for i in range(n_iteracoes):
        r1 = np.random.uniform(0, 1, (n_particulas, n_pgs+n_tensoes))
        r2 = np.random.uniform(0, 1, (n_particulas, n_pgs+n_tensoes))
        velocidades = w * velocidades + c1 * r1 * (pbest - posicoes) + c2 * r2 * (posicao_gbest - posicoes)
        posicoes += velocidades
        for j in range(n_particulas):
            posicoes[j, :n_tensoes] = np.clip(posicoes[j, :n_tensoes], Vmin, Vmax)
            posicoes[j, n_tensoes:] = np.clip(posicoes[j, n_tensoes:], Pgmin, 999999)
        fitness = np.array([MaxLoad(tipo, posicao, ybus) for posicao in posicoes])
        for j in range(n_particulas):
            if fitness[j] > valor_pbest[j]:
                valor_pbest[j] = fitness[j]
                pbest[j] = np.copy(posicoes[j])
        melhor_fitness_idx = np.argmax(fitness)
        if fitness[melhor_fitness_idx] > valor_gbest:
            valor_gbest = fitness[melhor_fitness_idx]
            posicao_gbest = np.copy(posicoes[melhor_fitness_idx])
        elif fitness[melhor_fitness_idx] == valor_gbest and performance(posicoes[melhor_fitness_idx]) < performance(posicao_gbest):
            posicao_gbest = np.copy(posicoes[melhor_fitness_idx])
        print(f'Iteração {i + 1} concluída! Carregamento máximo: {valor_gbest:.4f}. Posição: {posicao_gbest}')
    return valor_gbest, posicao_gbest

tipo = 'ATIVO/REATIVO'
n_pgs = 1
n_tensoes = 5

Pgmin = 0
Pgmax = 10
Vmin = 0.95
Vmax = 1.10

w = 0.8
c1 = 1.25
c2 = 1.25

n_particulas = 10
n_iteracoes = 20
n_repeticoes = 1

variaveis = ['V1', 'V2', 'V3', 'V6', 'V8', 'Pg2']

for i in range(n_repeticoes):
    t1 = time()
    carregamento_max, valor_vpg = pso(n_particulas, n_iteracoes, n_pgs, n_tensoes)
    t2 = time()

    print(f'Resultado da repetição {i+1}')
    for i in range(len(variaveis)):
        print(f'{variaveis[i]} = {valor_vpg[i]:.4f} pu')
    print(f'Carregamento máximo: {carregamento_max:.4f}')
    print(f'Tempo de execução: {t2-t1:.2f} s')



