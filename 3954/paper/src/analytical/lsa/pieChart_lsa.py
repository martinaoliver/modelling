import matplotlib.pyplot as plt
def pieChart_lsaf(valueCounts_dict,title):
    colors=['grey','grey','grey','peachpuff','coral','coral','coral','coral','coral','grey','coral','coral']
    labels = []
    sizes = []

    for x, y in valueCounts_dict.items():
        labels.append(x)
        sizes.append(y)

    plt.pie(sizes,colors=colors, labels=labels)
    plt.axis('equal')
    plt.title(title)
    plt.show()