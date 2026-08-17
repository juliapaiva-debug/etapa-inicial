# Preparo dos sistemas

## Usando o terminal do linux pela primeira vez

1. Abra o terminal do Linux e digite: 
```bash
pwd
```
Isso mostra o "endereço" completo da pasta onde você está agora (algo como /home/usuario/Documentos).

(aperte 'Enter' depois de cada comando)


Depois, digite:

```bash
ls
```
Isso mostra o que existe dentro dessa pasta (outros arquivos e pastas).

2. Se você não estiver na pasta certa, digite:

```bash
cd ..
```

Isso faz você subir um nível, saindo da pasta atual e indo para a pasta que a contém. Pode repetir esse comando quantas vezes precisar até chegar mais perto de onde quer.

Para entrar em uma pasta específica (por exemplo, uma chamada "projetos"), digite:
```bash
cd projetos
```
trocando "projetos" pelo nome real da pasta.

## Etapa 0 - Separação das cadeias

Antes de tudo, no caso da NS2B/NS3 com o ligante (6MO2), deve-se separar cada proteína nas cadeias A-C e B-D.

Para isso rodar o comando:

```bash
python split_ns2b_ns3.py 6MO2.pdb 6MO2_sep.pdb
```

## 📉 Etapa 1 - Separação das proteínas e verificação dos resíduos
A primeira etapa é separar as proteínas em um único arquivo para cada. Então, no pymol, selecionar "chains" e clicar nas cadeias que voce quer excluir. Em seguida, clicar em (sele) ➜ "A" ➜ "remove atoms". Restará somente a proteina com as cadeias de interesse. Salvar como denv_B_D. 

Conferir se há resíduos faltantes (sem densidade eletrônica) no seu sistema. Para isso, buscar por linhas tracejadas no seu .pdb aberto no Pymol, ou na aba "structure" no PDB. Caso haja, iremos modelar na proxima etapa. 

Outra parte importante, é a necessidade de verificar se os resíduos no .pdb estão com a numeração correta. Para isso, ir na aba "sequence" no PDB de sua estrutura. Nessa aba, havera o codigo Uniprot dessa respectiva proteina. Colocar esse codigo no Uniprot e ir ate a sessão "sequence", onde terá a sequência completa da proteina. Então, comparar a sequencia do seu .pdb com a sequencia de referencia da proteina e ver se a numeração esta compativel.

Dica: no pymol, ver a sequência dos 5 primeiros resíduos e dar cntrl + F no uniprot com esses 5 resíduos para encontra-los.

**OBS.:** Pode ocorrer que um desses 5 resíduos tenha um dos aminoácidos diferente; não se preocupe, o importante é que a maioria esteja igual ao do PDB aberto no pymol.

Caso a numeração não esteja batendo, será necessário reordenar os resíduos. Para isso, faremos um script no claude, nao se perde tempo com isso manualmente. 😎👍


## 📉 Etapa 2 - Modelagem comparativa
**Essa etapa requer a instalação do software Modeller**

Para modelar a região faltante, temos que modelar a sequência da nossa proteína contra o próprio modelo .pdb. 

Para obter a nossa sequência que será modelada, ir na página pdb da estrutura, > "Display" >"FASTA" > salvar essa sequência que será utilizada em breve.

**OBS.:** Caso haja sequencias adicionais em alguma das cadeias, como o linker de glicina, não é necessário copiá-lo.

**Parte 1 - Criação do arquivo .ali**

O arquivo .ali contém o seguinte formato:

```bash
>P1;denv
sequence:denv:::::::0.00: 0.00
GSHMLMADLELERAADVRWEEQAEISGS./
LEDGAYRIKQKGILGYSQIGAGVYKEGTFHTMWHVTRGAVLMHKGKRIEPSWADVKKDLISYGGGWKLEGEWKEGEEVQVLALEPGKNPRAVQTKPGLFKTNTGTIGAVSLDFSPGTSGSPIVDKKGKVVGLYGNGVVTRSGAYVSAIANTEKSIEDNPEIEDDIFRK*
```

Onde, na primeira linha deverá conter >P1;"nome qualquer"

Na segunda linha: "sequence:"mesmo nome":::::::0.00: 0.00"

A partir da terceira linha, deve-se conter a sequência que será modelada

**⚠️ O símbolo / representa uma quebra de cadeia, sempre deverá estar entre o ultimo residuo de uma cadeia e o primeiro da outra**

**⚠️ O símbolo . representa um ligante, sempre deverá estar na cadeia em que esta o ligante**

**⚠️ o símbolo * sempre será o último caracter da sequência**

**Parte 2 - Alinhamento entre sequência e o modelo**

O linhamento é realizado através do seguinte script:

```bash
from modeller import *

env = Environ()
env.io.hetatm = True       # ESSENCIAL: sem isso o Modeller ignora o ligante (HETATM)

aln = Alignment(env)

# --- Estrutura molde (template) ---
# troque 'template' pelo nome do seu PDB molde (sem extensão) e 'template.pdb' pelo arquivo real
mdl = Model(env, file='6MO2_A_C', model_segment=('FIRST:A', 'LAST:C')) #Colocar o nome correto das cadeias
aln.append_model(mdl, align_codes='templateAC', atom_files='6MO2_A_C.pdb')

# --- Sua sequência alvo ---
aln.append(file='denv.ali', align_codes='denv') #Colocar o nome o arquivo .ali e o nome dado para a sequência

# --- Alinhamento estrutura-sequência ---
aln.align2d(max_gap_length=50)

aln.write(file='denv-template.ali', alignment_format='PIR') #nomes dos outputs
aln.write(file='denv-template.pap', alignment_format='PAP')
```

Com isso, serão gerados os arquivos .pir e .pap de Alinhamento

**Parte 3 - Geração dos modelos**
A última parte no modeller é a geração dos modelos por homologia.

Para isso, usar o script:

```bash
from modeller import *
from modeller.automodel import *

env = Environ()
env.io.hetatm = True        # necessário para o ligante ser incluído no modelo final
# env.io.water = False       # opcional, se quiser excluir aguas

a = AutoModel(env, alnfile='denv-template.ali',
              knowns='templateAC', sequence='denv',  #Os nomes devem bater com o script anterior
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 100 #Numero de modelos a ser gerado
a.make()
```

**Parte 4 - Selecionando o melhor modelo**

O melhor modelo será aquele com o menor DOPE Score. Para descobrir de maneira automatizada, utilizar o comando:

```bash
grep denv.B9999 automodel.log | sort -n -k3
```

O melhor modelo será o primeiro da lista

## 📉 Etapa 3 - Refinando o modelo

Muitas vezes, o modelo gerado por modelagem comparativa não possui uma boa qualidade estrutural. Por isso, precisamos refiná-lo. 

A primeira etapa é determinar os estados de protonação dos residuos. Inicialmente, deve-se adicionar o .pdb modelado no PDB2PQR. Selecionar o pH 7.4, campo de força AMBER e saída CHARMM. Ao término do programa, baixar o arquivo .pqr.

Para verificar os estados de protonação, utilizar os comandos:

```bash
grep HSE arquivo.pqr
grep HSP arquivo.pqr
grep ASP arquivo.pqr
grep GLU arquivo.pqr
```
