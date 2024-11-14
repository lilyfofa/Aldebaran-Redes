from funcoes import Load
import numpy as np
import matplotlib.pyplot as plt

# def Load(alpha, opcao, slim=None, ybus=None, V=None, Pg=None, Qg=None, Ql=None, tipo=None):

n = 14

# FPO Ativo
print('FPO Ativo')
print('-='*50)
alpha = 0
step = 0.01

sucesso = True

ativo_y = []
ativo_x = []

while sucesso:
    valor = Load(alpha, True, True, Pg=np.array([0, 0.64194659, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64))
    sucesso = valor[0]
    if sucesso:
        print(f'Sucesso para alpha = {alpha:.4f}')
        ativo_x.append(alpha)
        ativo_y.append(valor[1])
        alpha += step
    else:
        print(f'Falha para alpha = {alpha:.4f}')

ativo_y = np.array(ativo_y)

# FPO reativo
print('FPO Reativo')
print('-='*50)
alpha = 0
step = 0.01

sucesso = True

reativo_y = []
reativo_x = []

while sucesso:
    valor = Load(alpha, True, True, V=np.array([1.1, 1.07704862, 0.98266485, 0, 0, 0.98156622, 0, 1.1, 0, 0, 0, 0, 0, 0], dtype=np.float64))
    sucesso = valor[0]
    if sucesso:
        print(f'Sucesso para alpha = {alpha:.4f}')
        reativo_x.append(alpha)
        reativo_y.append(valor[1])
        alpha += step
    else:
        print(f'Falha para alpha = {alpha:.4f}')

reativo_y = np.array(reativo_y)

# FPO Ativo-Reativo
print('FPO Ativo/Reativo')
print('-='*50)
alpha = 0
step = 0.01

sucesso = True

ativo_reativo_y = []
ativo_reativo_x = []

while sucesso:
    valor = Load(alpha, True, True, V=np.array([1.09475372, 1.1, 1.00100171, 0, 0, 0.95, 0, 1.00449198, 0, 0, 0, 0, 0, 0], dtype=np.float64),
                 Pg=np.array([0, 3.39757119, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64))
    sucesso = valor[0]
    if sucesso:
        print(f'Sucesso para alpha = {alpha:.4f}')
        ativo_reativo_x.append(alpha)
        ativo_reativo_y.append(valor[1])
        alpha += step
    else:
        print(f'Falha para alpha = {alpha:.4f}')

ativo_reativo_y = np.array(ativo_reativo_y)

print(len(ativo_x))
print(len(reativo_x))
print(len(ativo_reativo_x))
print(ativo_y.shape)
print(reativo_y.shape)
print(ativo_reativo_y.shape)

plt.figure(1)
for i in range(n):
    plt.plot(ativo_x, ativo_y[:, i], label=f'Barra {i+1}')
plt.title('Tensão nas barras versus carregamento (FPO ativo).')
plt.ylabel('Magnitude de tensão (p.u.)')
plt.xlabel('Carregamento')
plt.legend()
plt.grid()
plt.show()

plt.figure(2)
for i in range(n):
    plt.plot(reativo_x, reativo_y[:, i], label=f'Barra {i+1}')
plt.title('Tensão nas barras versus carregamento (FPO reativo).')
plt.ylabel('Magnitude de tensão (p.u.)')
plt.xlabel('Carregamento')
plt.legend()
plt.grid()
plt.show()

plt.figure(3)
for i in range(n):
    plt.plot(ativo_reativo_x, ativo_reativo_y[:, i], label=f'Barra {i+1}')
plt.title('Tensão nas barras versus carregamento (FPO ativo/reativo).')
plt.ylabel('Magnitude de tensão (p.u.)')
plt.xlabel('Carregamento')
plt.legend()
plt.grid()
plt.show()