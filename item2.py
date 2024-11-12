import numpy as np
from tabulate import tabulate

# Início
print('Nicolas melhor vice da empeltec')

# DCTE
BASE = 100
TEP = 0.01
erro = TEP/BASE

# DBAR
Pg = np.array([0, 0.40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64)
Qg = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64)
Pl = np.array([0, 0.217, 0.942, 0.478, 0.076, 0.112, 0, 0, 0.295, 0.09, 0.035, 0.061, 0.135, 0.149], dtype=np.float64)
Ql = np.array([0, 0.127, 0, -0.039, 0.016, 0.06, 0, 0, 0.166, 0.058, 0.018, 0.016, 0.058, 0.05], dtype=np.float64)
tipo = np.array(['SLACK', 'PV', 'PQ', 'PQ', 'PQ', 'PQ', 'PQ', 'PV', 'PQ', 'PQ', 'PQ', 'PQ', 'PQ', 'PQ'])

P = Pg - Pl
Q = Qg - Ql

V = np.array([1.02, 1.05, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0], dtype=np.float64)
θ = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64)
n = 14

pv, pq, slack = [], [], []
for i, tipo_barra in enumerate(tipo):
        if tipo_barra == 'PV':
            pv.append(i)
        elif tipo_barra == 'PQ':
            pq.append(i)
        elif tipo_barra == 'SLACK':
            slack.append(i)

# DLIN
dlin = np.array([
    (1, 2, 1.938, 5.917, 5.28),
    (1, 5, 5.403, 22.304, 4.92),
    (2, 3, 4.699, 19.797, 4.38),
    (2, 4, 5.811, 17.632, 3.74),
    (2, 5, 5.695, 17.388, 3.4),
    (3, 4, 6.701, 17.103, 3.56),
    (4, 5, 1.335, 4.211, 1.28),
    (4, 7, 0, 20.912, 0),
    (4, 9, 0, 55.618, 0),
    (5, 6, 0, 25.202, 0),
    (6, 11, 9.498, 19.89, 0),
    (6, 12, 12.291, 25.581, 0),
    (6, 13, 6.615, 13.027, 0),
    (7, 8, 0, 17.615, 0),
    (7, 9, 0, 11.001, 0),
    (9, 10, 3.181, 8.45, 0),
    (9, 14, 12.711, 27.038, 0),
    (10, 11, 8.205, 19.207, 0),
    (12, 13, 22.092, 19.988, 0),
    (13, 14, 17.093, 34.802, 0)], dtype=np.float64)

# Matrizes Y, G e B
ybus = np.zeros((n, n), dtype=complex)

for barra1, barra2, r, x, q in dlin:
    barra1 = int(barra1)
    barra2 = int(barra2)
    ybus[barra1-1][barra1-1] += 1/(r/100 + 1j * x/100) + 1j * q/(BASE*2)
    ybus[barra2-1][barra2-1] += 1/(r/100 + 1j * x/100) + 1j * q/(BASE*2)
    ybus[barra1-1][barra2-1] -= 1/(r/100 + 1j * x/100)
    ybus[barra2-1][barra1-1] -= 1/(r/100 + 1j * x/100)

# Compensador
Sh = 19
ybus[8][8] += Sh/BASE

# Display da matriz Ybus
lista_ybus = []
for i in range(0, n):
    linha = []
    for j in range(0, n):
        linha.append(complex(np.round(ybus[i][j], 4)))
    lista_ybus.append(linha)


print('-'*200)
print('Tabela 1 - YBUS do sistema-teste de 6 barras')
print(tabulate(lista_ybus, tablefmt='fancy_grid', numalign="center", stralign='center'))
print('-'*200)

# Identificação das variáveis
n_θ = 0
n_V = 0
index_θ = []
index_V = []

for i in range(1, n):
    if θ[i] == 0:
        n_θ += 1
        index_θ.append(i)
for i in range(0, n):
    if V[i] == 0:
        n_V += 1
        index_V.append(i)
        V[i] = 1

# Tabela das iterações, ajuste e fluxo
tabela_iteracoes = [['Iteração']]
tabela_ajuste = [['Iteração']]
tabela_fluxo = []
for index in index_θ:
    tabela_iteracoes[0].append(f'θ{index+1} [rad]')
    tabela_iteracoes[0].append(f'θ{index+1} [°]')
    tabela_ajuste[0].append(f'Δθ{index + 1} [rad]')
for index in index_V:
    tabela_iteracoes[0].append(f'V{index + 1} [pu]')
    tabela_ajuste[0].append(f'ΔV{index + 1} [pu]')

# Laço iterativo
iterations = 0
while True:
    # Equações ΔP
    ΔP = []
    index_ΔP = []
    for k in pv+pq:
        #if P[k] != 0:
            expr = P[k]
            for m in range(n):
                expr -= V[k]*V[m]*(np.real(ybus[k][m])*np.cos(θ[k]-θ[m]) + np.imag(ybus[k][m])*np.sin(θ[k]-θ[m]))
            ΔP.append(expr)
            index_ΔP.append(k)

    # Equações ΔQ
    ΔQ = []
    index_ΔQ = []
    for k in pq:
        #if Q[k] != 0:
            expr = Q[k]
            for m in range(n):
                expr -= V[k]*V[m]*(np.real(ybus[k][m])*np.sin(θ[k]-θ[m]) - np.imag(ybus[k][m])*np.cos(θ[k]-θ[m]))
            ΔQ.append(expr)
            index_ΔQ.append(k)

    # Vetor B
    B = np.concatenate((ΔP, ΔQ))
    linha_fluxo = list(B)
    linha_fluxo.insert(0, iterations)
    tabela_fluxo.append(linha_fluxo)

    # Critério de parada
    if max(abs(B)) < erro:
        linha = [iterations]
        for index in index_θ:
            linha.append(θ[index])
            linha.append(np.degrees(θ[index]))
        for index in index_V:
            linha.append(V[index])
        tabela_iteracoes.append(linha)
        break

    # Submatriz H = dΔP/dθ
    H = np.zeros((len(ΔP), n_θ))
    i = 0

    for k in index_ΔP:
        j = 0
        for m in index_θ:
            if k != m:
                H[i][j] = V[k]*V[m]*(np.real(ybus[k][m])*np.sin(θ[k]-θ[m]) - np.imag(ybus[k][m])*np.cos(θ[k]-θ[m]))
            else:
                expr = -V[k]**2*np.imag(ybus[k][k])
                for o in range(n):
                    expr -= V[k]*V[o]*(np.real(ybus[k][o])*np.sin(θ[k]-θ[o]) - np.imag(ybus[k][o])*np.cos(θ[k]-θ[o]))
                H[i][j] = expr
            j += 1
        i += 1

    # Submatriz N = dΔP/dV
    N = np.zeros((len(ΔP), n_V))
    i = 0

    for k in index_ΔP:
        j = 0
        for m in index_V:
            if k != m:
                N[i][j] = V[k]*(np.real(ybus[k][m])*np.cos(θ[k]-θ[m]) + np.imag(ybus[k][m])*np.sin(θ[k]-θ[m]))
            else:
                expr = V[k]*np.real(ybus[k][k])
                for o in range(n):
                    expr += V[o]*(np.real(ybus[k][o])*np.cos(θ[k]-θ[o]) + np.imag(ybus[k][o])*np.sin(θ[k]-θ[o]))
                N[i][j] = expr
            j += 1
        i += 1

    # Submatriz M = dΔQ/dθ
    M = np.zeros((len(ΔQ), n_θ))
    i = 0

    for k in index_ΔQ:
        j = 0
        for m in index_θ:
            if k != m:
                M[i][j] = -V[k]*V[m]*(np.real(ybus[k][m])*np.cos(θ[k]-θ[m]) + np.imag(ybus[k][m])*np.sin(θ[k]-θ[m]))
            else:
                expr = -V[k]**2*np.real(ybus[k][k])
                for o in range(n):
                    expr += V[k]*V[o]*(np.real(ybus[k][o])*np.cos(θ[k]-θ[o]) + np.imag(ybus[k][o])*np.sin(θ[k]-θ[o]))
                M[i][j] = expr
            j += 1
        i += 1

    # Submatriz L = dΔQ/dV
    L = np.zeros((len(ΔQ), n_V))
    i = 0

    for k in index_ΔQ:
        j = 0
        for m in index_V:
            if k != m:
                L[i][j] = V[k]*(np.real(ybus[k][m])*np.sin(θ[k]-θ[m]) - np.imag(ybus[k][m])*np.cos(θ[k]-θ[m]))
            else:
                expr = -V[k]*np.imag(ybus[k][k])
                for o in range(n):
                    expr += V[o]*(np.real(ybus[k][o])*np.sin(θ[k]-θ[o]) - np.imag(ybus[k][o])*np.cos(θ[k]-θ[o]))
                L[i][j] = expr
            j += 1
        i += 1

    # Jacobiano
    J = np.block([[H, N], [M, L]])

    # Ajuste
    Δx = np.linalg.solve(J, B)
    linha_ajuste = list(Δx)
    linha_ajuste.insert(0, iterations)

    tabela_ajuste.append(linha_ajuste)

    # Linha da tabela das iterações
    linha = [iterations]

    contador = 0

    for index in index_θ:
        linha.append(θ[index])
        linha.append(np.degrees(θ[index]))
        θ[index] += Δx[contador]
        contador += 1

    for index in index_V:
        linha.append(V[index])
        V[index] += Δx[contador]
        contador += 1

    tabela_iteracoes.append(linha)

    iterations += 1

tabela_J = []

for index in index_ΔP:
    linha = []
    for index2 in index_θ:
        linha.append(f'dP{index+1}/dθ{index2+1}')
    for index2 in index_V:
        linha.append(f'dP{index+1}/dV{index2+1}')
    tabela_J.append(linha)

for index in index_ΔQ:
    linha = []
    for index2 in index_θ:
        linha.append(f'dQ{index+1}/dθ{index2+1}')
    for index2 in index_V:
        linha.append(f'dQ{index+1}/dV{index2+1}')
    tabela_J.append(linha)

print("Tabela 2 - Matriz jacobiana simbólica do sistema-teste.")
print(tabulate(tabela_J, headers='firstrow', tablefmt='fancy_grid', numalign="center", floatfmt=".4f"))
print("-"*200)

print("Tabela 3 - Processo iterativo de resolução do problema.")
print(tabulate(tabela_iteracoes, headers='firstrow', tablefmt='fancy_grid', numalign="center", floatfmt=".4f"))
print("-"*200)

header_fluxo = ['Iteração']
for index in index_ΔP:
    header_fluxo.append(f'ΔP{index+1} [pu]')
for index in index_ΔQ:
    header_fluxo.append(f'ΔQ{index+1} [pu]')
tabela_fluxo.insert(0, header_fluxo)

print("Tabela 4 - Evolução dos fluxos de potência ao longo do processo iterativo.")
print(tabulate(tabela_fluxo, headers='firstrow', tablefmt='fancy_grid', numalign="center", floatfmt=".4f"))
print("-"*200)

print("Tabela 5 - Ajuste das variáveis ao longo do processo iterativo.")
print(tabulate(tabela_ajuste, headers='firstrow', tablefmt='fancy_grid', numalign="center", floatfmt=".4f"))
print("-"*200)

# Dados de barra após resolução do sistema
tabela_dbar = [['Barra', 'V [pu]', 'θ [rad]', 'θ [°]', 'P [pu]', 'Q [pu]']]

# Encontrando Pk e Qk

for k in range(0, n):
    if P[k] == 0:
        expressao = 0
        for m in range(0, n):
            expressao += V[k]*V[m]*(np.real(ybus[k][m])*np.cos(θ[k]-θ[m]) + np.imag(ybus[k][m])*np.sin(θ[k]-θ[m]))
        P[k] = expressao
    if Q[k] == 0:
        expressao = 0
        for m in range(0, n):
            expressao += V[k]*V[m]*(np.real(ybus[k][m])*np.sin(θ[k]-θ[m]) - np.imag(ybus[k][m])*np.cos(θ[k]-θ[m]))
        Q[k] = expressao
    tabela_dbar.append([k+1, V[k], θ[k], np.degrees(θ[k]), P[k], Q[k]])

print("Tabela 6 - Resumo dos dados de barra após resolução do sistema.")
print(tabulate(tabela_dbar, headers='firstrow', tablefmt='fancy_grid', numalign="center", floatfmt=".4f"))
print("-"*200)

tabela_dlin = [['Linha', 'R [pu]', 'X [pu]', 'Bsh [pu]', 'Pkm [pu]', 'Pmk [pu]', 'Qkm [pu]', 'Skm [pu]', 'Perdas [pu]']]

soma = 0
for barra1, barra2, r, x, q in dlin:
    k = int(barra1) - 1
    m = int(barra2) - 1
    rkm = r/100
    xkm = x/100
    gkm = rkm/(np.square(rkm)+np.square(xkm))
    bkm = -xkm/(np.square(xkm)+np.square(rkm))
    bsh = q/(2*BASE)
    pkm = np.square(V[k])*gkm - V[k]*V[m]*gkm*np.cos(θ[k]-θ[m]) - V[k]*V[m]*bkm*np.sin(θ[k]-θ[m])
    pmk = np.square(V[m])*gkm - V[k]*V[m]*gkm*np.cos(θ[k]-θ[m]) + V[k]*V[m]*bkm*np.sin(θ[k]-θ[m])
    perdas = pkm + pmk
    qkm = -np.square(V[k])*(bkm+bsh) + V[k]*V[m]*bkm*np.cos(θ[k]-θ[m]) - V[k]*V[m]*gkm*np.sin(θ[k]-θ[m])
    skm = np.sqrt(pkm**2 + qkm**2)
    tabela_dlin.append([f'{barra1}-{barra2}', r/100, x/100, bsh, pkm, pmk, qkm, skm, perdas])
    soma += perdas

print("Tabela 7 - Resumo dos dados de linha após resolução do sistema.")
print(tabulate(tabela_dlin, headers='firstrow', tablefmt='fancy_grid', numalign="center", floatfmt=".4f"))
print(f'Perdas totais: {soma:.4f} pu')
print("-"*200)

print("(c) 2024 - EmpelTec Jr. - Todos os direitos reservados")