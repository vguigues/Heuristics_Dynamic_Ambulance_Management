import pandas as pd
import matplotlib.pyplot as plt

data = {
    'A': [10, 12, 14, 16, 18],
    'avg enum call': [454.7, 521.0, 672.4, 672.9, 699.1],
    'std enum call': [423.7, 612.3, 859.8, 890.4, 901.1],
    'avg enum amb': [463.0, 536.7, 679.9, 698.0, 723.3],
    'std enum amb': [616.3, 1074.6, 1101.0, 1719.1, 1679.8],
    'avg call': [1791.6, 1918.4, 2328.1, 991.7, 1042.6],
    'std call': [398.5, 263.9, 382.2, 300.9, 286.0],
    'avg amb': [1720.5, 1823.1, 2277.4, 905.9, 963.0],
    'std amb': [363.5, 316.9, 370.7, 293.7, 339.5]
}

df = pd.DataFrame(data)

fig, ax = plt.subplots(figsize=(10, 6))

x = df['A']

for i in range(1, 9, 2):
    y = df.iloc[:, i]
    y_err = df.iloc[:, i+1]
    ax.plot(x, y, label=df.columns[i])
    # ax.fill_between(x, y - y_err, y + y_err, alpha=0.3)

ax.legend()
ax.set_xlabel('A')
ax.set_ylabel('Valor')
ax.set_title('Média e desvio padrão para cada par de colunas')
ax.set_ylim(bottom=0)

plt.show()