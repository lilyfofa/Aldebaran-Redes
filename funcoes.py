import numpy as np


def Load(alpha, opcao, V=None, Pg=None, Qg=None, Ql=None, tipo=None):
    # Definição dos valores padrões
    if V is None:
        V = np.array([1.02, 1.05, 1.00, 0, 0, 1.00, 0, 1.00, 0, 0, 0, 0, 0, 0], dtype=np.float64)
    if Pg is None:
        Pg = np.array([0, 0.40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64)
    if Qg is None:
        Qg = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64)
    if Ql is None:
        Ql = np.array([0, 0.127, 0.19, -0.039, 0.016, 0.075, 0, 0, 0.166, 0.058, 0.018, 0.016, 0.058, 0.05],
            dtype=np.float64)
    if tipo is None:
        tipo = np.array(['SLACK', 'PV', 'PV', 'PQ', 'PQ', 'PV', 'PQ', 'PV', 'PQ', 'PQ', 'PQ', 'PQ', 'PQ', 'PQ'])
    # DCTE
    BASE = 100
    TEP = 0.01
    erro = TEP / BASE

    # DBAR
    Pl = np.array([0, 0.217, 0.942, 0.478, 0.076, 0.112, 0, 0, 0.295, 0.09, 0.035, 0.061, 0.135, 0.149],
                  dtype=np.float64)

    P = Pg * (1 + alpha) - Pl * (1 + alpha)
    Q = Qg - Ql * (1 + alpha)

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
        ybus[barra1 - 1][barra1 - 1] += 1 / (r / 100 + 1j * x / 100) + 1j * q / (BASE * 2)
        ybus[barra2 - 1][barra2 - 1] += 1 / (r / 100 + 1j * x / 100) + 1j * q / (BASE * 2)
        ybus[barra1 - 1][barra2 - 1] -= 1 / (r / 100 + 1j * x / 100)
        ybus[barra2 - 1][barra1 - 1] -= 1 / (r / 100 + 1j * x / 100)

    # Compensador
    Sh = 19
    ybus[8][8] += Sh / BASE

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

    # Laço iterativo
    iterations = 0
    while True:
        # Equações ΔP
        ΔP = []
        index_ΔP = []
        for k in pv + pq:
            # if P[k] != 0:
            expr = P[k]
            for m in range(n):
                expr -= V[k] * V[m] * (
                        np.real(ybus[k][m]) * np.cos(θ[k] - θ[m]) + np.imag(ybus[k][m]) * np.sin(θ[k] - θ[m]))
            ΔP.append(expr)
            index_ΔP.append(k)

        # Equações ΔQ
        ΔQ = []
        index_ΔQ = []
        for k in pq:
            # if Q[k] != 0:
            expr = Q[k]
            for m in range(n):
                expr -= V[k] * V[m] * (
                        np.real(ybus[k][m]) * np.sin(θ[k] - θ[m]) - np.imag(ybus[k][m]) * np.cos(θ[k] - θ[m]))
            ΔQ.append(expr)
            index_ΔQ.append(k)

        # Vetor B
        B = np.concatenate((ΔP, ΔQ))

        # Critério de parada
        if max(abs(B)) < erro:
            break

        # Submatriz H = dΔP/dθ
        H = np.zeros((len(ΔP), n_θ))
        i = 0

        for k in index_ΔP:
            j = 0
            for m in index_θ:
                if k != m:
                    H[i][j] = V[k] * V[m] * (
                            np.real(ybus[k][m]) * np.sin(θ[k] - θ[m]) - np.imag(ybus[k][m]) * np.cos(θ[k] - θ[m]))
                else:
                    expr = -V[k] ** 2 * np.imag(ybus[k][k])
                    for o in range(n):
                        expr -= V[k] * V[o] * (np.real(ybus[k][o]) * np.sin(θ[k] - θ[o]) - np.imag(ybus[k][o]) * np.cos(
                            θ[k] - θ[o]))
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
                    N[i][j] = V[k] * (
                            np.real(ybus[k][m]) * np.cos(θ[k] - θ[m]) + np.imag(ybus[k][m]) * np.sin(θ[k] - θ[m]))
                else:
                    expr = V[k] * np.real(ybus[k][k])
                    for o in range(n):
                        expr += V[o] * (np.real(ybus[k][o]) * np.cos(θ[k] - θ[o]) + np.imag(ybus[k][o]) * np.sin(
                            θ[k] - θ[o]))
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
                    M[i][j] = -V[k] * V[m] * (
                            np.real(ybus[k][m]) * np.cos(θ[k] - θ[m]) + np.imag(ybus[k][m]) * np.sin(θ[k] - θ[m]))
                else:
                    expr = -V[k] ** 2 * np.real(ybus[k][k])
                    for o in range(n):
                        expr += V[k] * V[o] * (np.real(ybus[k][o]) * np.cos(θ[k] - θ[o]) + np.imag(ybus[k][o]) * np.sin(
                            θ[k] - θ[o]))
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
                    L[i][j] = V[k] * (
                            np.real(ybus[k][m]) * np.sin(θ[k] - θ[m]) - np.imag(ybus[k][m]) * np.cos(θ[k] - θ[m]))
                else:
                    expr = -V[k] * np.imag(ybus[k][k])
                    for o in range(n):
                        expr += V[o] * (np.real(ybus[k][o]) * np.sin(θ[k] - θ[o]) - np.imag(ybus[k][o]) * np.cos(
                            θ[k] - θ[o]))
                    L[i][j] = expr
                j += 1
            i += 1

        # Jacobiano
        J = np.block([[H, N], [M, L]])

        # Ajuste
        Δx = np.linalg.solve(J, B)

        contador = 0

        for index in index_θ:
            θ[index] += Δx[contador]
            contador += 1

        for index in index_V:
            V[index] += Δx[contador]
            contador += 1

        iterations += 1

        # Critério de parada
        if iterations > 50:
            return [False]
    # Encontrando Pk e Qk
    for k in range(0, n):
        if P[k] == 0:
            expressao = 0
            for m in range(0, n):
                expressao += V[k] * V[m] * (
                        np.real(ybus[k][m]) * np.cos(θ[k] - θ[m]) + np.imag(ybus[k][m]) * np.sin(θ[k] - θ[m]))
            P[k] = expressao
        if Q[k] == 0:
            expressao = 0
            for m in range(0, n):
                expressao += V[k] * V[m] * (
                        np.real(ybus[k][m]) * np.sin(θ[k] - θ[m]) - np.imag(ybus[k][m]) * np.cos(θ[k] - θ[m]))
            Q[k] = expressao

    # Verificando violações de limites
    if opcao:
        Qn = (1 / BASE) * np.array([-9999, -40, 0, 0, 0, -6, 0, -6, 0, 0, 0, 0, 0, 0], dtype=np.float64)
        Qm = (1 / BASE) * np.array([500, 80, 70, 0, 0, 40, 0, 40, 0, 0, 0, 0, 0, 0], dtype=np.float64)
        violacao = False
        V_novo = np.array([1.02, 1.05, 1.00, 0, 0, 1.00, 0, 1.00, 0, 0, 0, 0, 0, 0], dtype=np.float64)
        Qg_novo = Qg
        Ql_novo = Ql
        tipo_novo = tipo
        for k in range(len(Q)):
            if (Q[k] > Qm[k] or Q[k] < Qn[k]) and tipo[k] == 'PV':
                if Q[k] > Qm[k]:
                    Qg_novo[k] = Qm[k]
                else:
                    Qg_novo[k] = Qn[k]
                V_novo[k] = 0
                Ql_novo[k] = 0
                tipo_novo[k] = 'PQ'
                violacao = True
        if violacao:
            return Load(alpha, False, V_novo, Pg, Qg_novo, Ql_novo, tipo_novo)
        else:
            return True, V
    else:
        return True, V


def MaxLoad(tipo, posicao):
    if tipo == 'ATIVO':
        a = 0
        b = 5
        while b - a > 0.0001:
            c = 0.5 * (a + b)
            valor_c = Load(c, True, Pg=np.array([0, posicao[0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64))
            # print(c, valor_c)
            if valor_c[0]:
                a = c
            else:
                b = c
        # print(a)
        return a

