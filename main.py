"""
22/02/2022
Projeto 1 - Teoria e Aplicação de Grafos - 2021.2
Universidade de Brasília - Instituto de Ciências Exatas - Departamento de Ciência da Computação (CiC)
Professor: Díbio Leandro Borges

O objetivo desse programa é aplicar o algoritimo de Bron-Kerbosh para achar cliques maximais em um grafo e encontrar o coeficiente de aglomeração. 
O grafo que foi utilizado foi disponibilizado no artigo “David Lusseau et al., The bottelenose dolphin community of Doubful Sound features a large proportion of long-lasting associations, Journal of Behavioral Ecology and Sociobiology 54:4, 396--405 (2003).”. 
Ele consiste em uma rede social de relações duradouras é identificada em uma comunidade de 62 golfinhos e apresentada como um grafo (não direcionado) para estudos. 
O programa lê o arquivo (soc-dolphins.mtx), monta com esses dados um grafo não direcionado, sem pesos, e então realiza o seguinte:
(1) Implementa duas formas do algoritmo Bron-Kerbosch: uma com pivotamento, outra sem pivotamento;
(2) Encontra e imprime na tela (duas vezes, uma para cada implementação do Bron-Kerbosch) todos os cliques maximais (indicando o número de vértices e quais);
(3) O Coeficiente médio de Aglomeração do Grafo.
"""

import random  

# Função de leitura e montagem da lista de adjacências
def read_adjacences_from_file(file_name):
    with open(file_name, 'r') as f: # Leitura do arquivo
        adjacences = [[] for i in range(62)] # Inicializa lista de adjacências

        # Lê as linhas e coleta apenas as informações relevantes, aplicando as conversões necessárias:
        # - Ignora o cabeçalho e a linha com número de vértices e tamanho;
        # - Converte cada linha para uma tupla que representa uma aresta;
        # - Subtrai 1 de todos os vértices, para encaixar na indexação do objeto do tipo lista (será somado novamente no final, antes de printar). 
        edges = [tuple(map(lambda x: int(x) - 1, i.split())) for i in f.readlines()[36:]]

        # Define a lista de adjacências baseando-se nas tuplas geradas (conectando cada vértice ao seu vizinho)
        for i in edges:
            adjacences[i[0]].append(i[1])
            adjacences[i[1]].append(i[0])

        return adjacences # Retorna a lista de adjacências



# Função de Bron-Kerbosh com pivoteamento
def bron_kerbosh_with_pivoting(adjacences):
    def internal_bron_kerbosh_with_pivoting(P, R, X): # Função interna para recursão
        if not P and not X: # Se P e X forem vazias...
            yield R # então R é um clique maximal.
                
        U = P | set(X) # União de P com X

        # Se a união de P com X for vazia, seleciona P, se não, seleciona um pivô aleatório dessa união...
        # ...e depois faz um loop chamando Bron-Kerbosh recursivamente para todos os vértices restantes que foram selecionados.
        for v in (P if not U else P - set(adjacences[random.choice(list(U))])):

            # Chama recursivamente a função, passando a lista de adjacência:
            # P é a intersecção do P anterior com os adjacentes do vértice atual;
            # R é a união do R anterior com o vértice atual;
            # X é a intersecção do X anterior com os adjacentes do vértice atual.
            yield from internal_bron_kerbosh_with_pivoting(P & set(adjacences[v]), R | {v}, X & set(adjacences[v]))

            P -= {v} # Retira o vértice atual de P
            X |= {v} # Adiciona o vértice atual ao X
    
    # Prepara as entradas iniciais
    init_P = set(range(len(adjacences))) # P = Vértices
    cliques = internal_bron_kerbosh_with_pivoting(init_P, set(), set()) # Chamada inicial da função recursiva
    return [{i + 1 for i in x} for x in cliques] # Adiciona novamente o 1 que foi retirado para a montagem da lista e retorna



# Função de Bron-Kerbosh sem pivoteamento
def bron_kerbosh_without_pivoting(adjacences):
    def internal_bron_kerbosh_without_pivoting(P, R, X): # Função interna para recursão
        if not P and not X: # Se P e X forem vazias...
            yield R # então R é um clique maximal.

        while P: # Seleciona um por um cada vértice até P ser vazio
            v = P.pop() # Retira o vértice atual de P

            # Chama recursivamente a função, passando a lista de adjacência:
            # P é a intersecção do P anterior com os adjacentes do vértice atual;
            # R é a união do R anterior com o vértice atual;
            # X é a intersecção do X anterior com os adjacentes do vértice atual.
            yield from internal_bron_kerbosh_without_pivoting(P & set(adjacences[v]), R | {v}, X & set(adjacences[v]))

            X |= {v} # Adiciona o vértice atual ao X
    
    # Prepara as entradas iniciais
    init_P = set(range(len(adjacences))) # P = Vértices
    cliques = internal_bron_kerbosh_without_pivoting(init_P, set(), set()) # Chamada inicial da função recursiva
    return [{i+1 for i in x} for x in cliques] # Adiciona novamente o 1 que foi retirado para a montagem da lista e retorna



# Função para calcular o coeficiente de aglomeração pela lista de adjacências
def clustering_coefficient(adjacences):
    length = len(adjacences) # número de vértices

    coeficients = 0 # Inicializa soma dos coeficientes locais
    for v in set(range(length)): # Loopar por cada vértice para calcular coeficientes
        adjacents = adjacences[v] # Pega os vértices adjacentes ao vértice atual
        num_of_adjacents = len(adjacents) # Número de adjacentes do atual

        temp = 0 # Inicia contador
        if num_of_adjacents > 1: # Caso exista mais de 1 adjacente, o denominador não será zero
            combinations = ((x, y) for y in adjacents for x in adjacents) # Todas as combinações possíveis entre os adjacentes
            numerator = sum(map(lambda i: 1 if i[1] in adjacences[i[0]] else 0, combinations)) # Calcula quantos são, de fato, vizinhos entre si (numerador)
            denominator = num_of_adjacents * (num_of_adjacents - 1) # Número de conexões possíveis entre os vizinhos (denominador)
            temp = numerator / denominator # Calcula proporção (coeficiente local em si)
        coeficients += temp # Adiciona coeficiente local ao total

    return coeficients / length # Retorna proporção (coeficiente global em si)



if __name__ == "__main__": # Inicia o projeto
    adjacences = read_adjacences_from_file('soc-dolphins.mtx') # Inicialização da lista de adjacências
    
    # Chama a função de Bron-Kerbosh com pivoteamento e printa cada clique em uma linha, separadamente
    print("Cliques maximais - Bron-Kerbosch com pivoteamento:")
    print(*bron_kerbosh_with_pivoting(adjacences), sep='\n')
    print()
    # Chama a função de Bron-Kerbosh sem pivoteamento e printa cada clique em uma linha, separadamente
    print("Cliques maximais - Bron-Kerbosch sem pivoteamento:")
    print(*bron_kerbosh_without_pivoting(adjacences), sep='\n')
    print()
    # Chama a função de calcular o coeficiente do grafo baseado na lista de adjacências
    print("Coeficiente de aglomeração:")
    print(clustering_coefficient(adjacences))