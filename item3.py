from funcoes import Load
import numpy as np
import matplotlib.pyplot as plt

n = 14

step = 0.01

alpha = 0

sucesso = True

y = []
x = []

while sucesso:
    resultado = Load(alpha, False)
    sucesso = resultado[0]
    if sucesso:
        print(f'Sucesso para alpha = {alpha:.4f}')
        y.append(resultado[1])
        x.append(alpha)
        alpha += step

matriz_Y = np.array(y)

plt.figure(1)
for i in range(n):
    plt.plot(x, matriz_Y[:, i], label=f'Barra {i+1}')
plt.title('Tensão nas barras versus carregamento sem considerar limites de Qg.')
plt.ylabel('Magnitude de tensão (p.u.)')
plt.xlabel('Carregamento')
plt.legend()
plt.grid()
plt.show()

