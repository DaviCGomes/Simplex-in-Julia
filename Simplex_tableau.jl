using LinearAlgebra

#Ler arquivo de texto
type, c, A, b, sinais = read_file("trabalho.txt")

#Verifica se algum elemento em 'b' é negativo
if !all(b .≥ 0)
    for i in eachindex(b)
        if b[i] ≥ 0
            continue
        end
        
        #Inverte o sinal de 'b' e de 'A'
        b[i] = -b[i]
        A[i] = -A[i]

        #Troca a desigualdade
        if sinais[i] == "<="
            sinais[i] == ">="
        elseif sinais[i] == ">="
            sinais[i] == "<="
        end
    end
end

if(type == "min")
    c = -c
end

#Pega o tamanho da matriz 'A'
m, n = size(A)

folga = zeros(m)
M = 0

for i in eachindex(sinais)
    if sinais[i] == "<="
        #Adiciona folga positiva
        if type == "max"
            aux = zeros(m)
            aux[i] = 1
            folga = [folga aux]
            M = [M 0]
        #Adiciona folga negativa, positiva e penalidade
        elseif type == "min"
            aux = zeros(m,2)
            aux[i,1] = -1
            aux[i,2] = 1
            folga = [folga aux]
            M = [M 0 1]
        end
    elseif  sinais[i] == ">="
        #Adiciona folga negativa, positiva e penalidade
        if type == "max"
            aux = zeros(m,2)
            aux[i,1] = -1
            aux[i,2] = 1
            folga = [folga aux]
            M = [M 0 1]
        #Adiciona folga positiva
        elseif type == "min"
            aux = zeros(m)
            aux[i] = 1
            folga = [folga aux]
            M = [M 0]
        end
    #Adiciona folga positiva e penalidade
    elseif  sinais[i] == "="
        aux = zeros(m)
        aux[i] = 1
        folga = [folga aux]
        if type == "max"
            M = [M 1]
        elseif type == "min"
            M = [M 1]
        end
    end
end

folga = folga[:, 2:end]
M = M[2:end]

#Adicionar as variáveis de folga em 'A' e 'c'
#Concatena a matriz 'A' com 'N'
A = [A folga]

#Atribui as variáveis de folga no 'c'
#folga = zeros(n-m, 1)
c = [c; M]
c = vec(c)

IB = create_base(A)
A[:,IB]
#with_logger(ConsoleLogger(stderr, Logging.Debug)) do
    z, c_final, A_final, b_final, IB, status, inter, Δt = simplex_tableau(c, A, b, type)    
#end

m, n = size(A)
tableau = [
    ' ' collect(1:n)' ' ';
    'z' c_final' z;
    IB A_final b_final
]

dual = A[:, IB]' \ c[IB]

bmax, bmin = sensitivity(A_final, b_final)

Δb = ["Increase" "Decrease";
    bmax bmin]


function simplex_tableau(c, A, b, type; max_inter = 100, max_time = 10.0, ϵ = sqrt(eps()))
    #Pega o tamano de 'A'
    m, n = size(A)
    
    if type == "max"
        c = -c
    end

    #Define a base
    IB = collect(n-m+1:n)
    @debug("", IB)

    #Cria a tableau
    tableau = [
        c'   0.0
        A   b
    ]
    @debug("", tableau)

    inter = 0
    start_time = time()
    Δt = 0.0

    status = :unknown
    solved = false
    tired = inter ≥ max_inter > 0 || Δt ≥ max_time > 0

    while !(solved || tired)
        @debug("Interação $inter")

        #Verifica se não há valores negativos em 'c'
        if all(c .≥ -ϵ)
            solved = true
            continue
        end

        #Escolha do elemento que entrará na base
        j = argmin(c)
        @debug("", j)

        #Pega todos os valores da coluna escolhida
        d = A[:, j]

        #Impede que entre valores negativos
        for i = 1:m
            if d[i] < -ϵ
            d[i] = 0.0
            end
        end

        #Divide a 'b' com a coluna escolhida e pega o menor resultado
        k = argmin(b ./ d)
        @debug("", k)

        #Atualiza a base
        IB[k] = j
        @debug("", IB)

        #Tranforma o elemento (k,j) do tableau em 1
        tableau[k+1, :] /= tableau[k+1, j]
        #Zera os elementos (i,j) do tableau, para k ≠ j
        for i = 1:m+1
            if i == k+1
                continue
            end
            tableau[i, :] -= tableau[k+1,:] * tableau[i,j]
        end

        @debug("", tableau)

        #Atualizar o 'c', 'A' e 'b' com os novos valores
        c = tableau[1, 1:end-1]
        @debug("", c)

        A = tableau[2:end, 1:end-1]
        @debug("", A)

        b = tableau[2:end, end]
        @debug("", b)

        inter += 1
        Δt = time() - start_time
        tired = inter ≥ max_inter > 0 || Δt ≥ max_time > 0
    end

    if solved
        status = :solved
    elseif tired
        if inter ≥ max_inter > 0
            status = :max_inter
        elseif Δt ≥ max_time > 0
            status = :max_time
        end
    end

    z = tableau[1, end]

    #A = A[sortperm(IB), :]
    #b = b[sortperm(IB)]
    #IB = sort(IB)

    if type == "min"
        z = -z
        A[:, m+1:end] = -A[:, m+1:end]
    end

    return z, c, A, b, IB, status, inter, Δt
    
end

function create_base(A)
    m,n = size(A)
    IB = zeros(m)
    j = 1

    for i = m+1:n
        o = count(x -> x == 0, A[:, i])
        l = count(x -> x == 1, A[:, i])
        if o+l == m
        IB[j] = i
        j += 1
        end
    end

    if !all(IB .> 0)
        for i = 1: count(x -> x == 0, IB)
            IB[j] = i
            j += 1
        end
    end

    return Int.(IB)
end

function sensitivity(A, b)
    m,n = size(A)
    
    Δb = (-b) ./ A[:,m+1:end]
    m,n = size(Δb)

    bmax = zeros(m)
    bmin = zeros(m)

    for i=1:m
        if all(Δb[:,i] .≥ 0)
            bmax[i] = findmin(Δb[findall(x-> x ≥ 0, Δb[:,i]),i])[1]
            bmin[i] = 0
        elseif all(Δb[:,i] .≤ 0)
            bmax[i] = 0
            bmin[i] = findmin(-Δb[findall(x-> x ≤ 0, Δb[:,i]),i])[1]
        else
            bmax[i] = findmin(Δb[findall(x-> x ≥ 0, Δb[:,i]),i])[1]
            bmin[i] = findmin(-Δb[findall(x-> x ≤ 0, Δb[:,i]),i])[1]
        end
    end

    return bmax, bmin
end

function read_file(file_name)
    #Lê o arquivo de texto
    lines = readlines(file_name)

    #Pega a primeira linha onde terá as informações do problema
    f = lines[1]

    #Pega as linhas com as restrições
    st = String[]
    st = fill("", size(lines)[1]-3)

    for i=3: size(lines)[1]-1
        st[i-2] = lines[i]
    end

    #Separa cada elemento de 'f'
    aux = split(f)

    #Pega o tipo do problema
    type = aux[1]

    #Cria o vertor 'c'
    c = get_elements(aux, size(aux)[1], 2)

    #Pega o número de 'x's
    j = size(c)[1]

    #Pega o número de restrinções
    m = size(st)[1]

    #Incializa as matrizes 'A', folga, 'b' e pega o sinal da inequação em cada linha
    A = zeros(m, j)
    sinais = fill("", m)
    b = zeros(m)

    for i = 1 : m
        aux = split(st[i])
        A[i,:] = get_elements(aux, size(aux)[1]-2, 1)
        b[i] = tryparse(Int, aux[end])
        sinais[i] = aux[end-1]
    end

    return type, c, A, b, sinais
end

#Recebe um vetor de entrada, o tamanho, e o inicio
#Retorna um vetor
function get_elements(str, s, start)
    #Incializa o vetor a ser retornado
    vec = zeros((s + 1) ÷ 3)

    #Variáveis para verificar se o número é negativo, a quantidade de 'x's e pegar os valores
    isnegative = false
    j = 1
    num = 0

    #Analiza cada elemento para a tipagem certa
    for i = start : s
        c = str[i]
        if c =="-"
            isnegative = true
        elseif c == "+"
            isnegative = false
        elseif c[1] == 'x'
            j = c[2] - '0'
            vec[j] = num
        else
            num = tryparse(Int, c)
        end
    end

    return vec
end

function export_result()
    
end