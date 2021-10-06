# Gibbs-sampler
A motif is defined as a sequence of nucleotides or proteins that has some specific biological function or structure. For example, the transcription factor binding sites (TFBS), as a kind of DNA motif, play an important role in regulating gene expression. Up to now, A large number of algorithms for finding DNA motifs have been developed. Gibbs sampler is one of the Monte Carlo algorithms for de novo motif discovery. De novo 
motif discovery refers to the process of finding a motif in a sequence but we do not know what the motif is like 
and where it is located in the sequence. 
### The main algorithm of Gibbs sampler is consist of following steps. 
#### 1. Initialization

- The first is random initialization in which a random set of values is assigned and select motifs from these points. So our first set of **_N_** (total number of the input sequences) "motifs" is essentially a random set of sequences of length **_L_** and is not expected to have any pattern. 

|  | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | C | A | A | C | G | G | T |
| 1 | G | A | T | T | C | C | G |
| 2 | C | G | C | G | G | A | G |
| 3 | A | T | C | C | T | C | C |
| 4 | T | T | A | C | G | C | T |
| 5 | G | T | A | A | A | G | G |
| 6 | C | A | T | G | A | C | A |
| 7 | T | C | A | T | G | A | C |
| 8 | T | A | A | C | A | C | A |
| 9 | A | C | A | C | C | C | A |
| 10 | A | A | G | G | C | C | A |
| 11 | A | A | A | A | G | C | G |
| 12 | G | A | G | A | T | C | T |
| 13 | C | G | C | G | A | G | G |
| 14 | G | G | G | C | A | C | T |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 98 | A | C | C | C | A | C | G |
| 99 | T | T | G | C | A | G | C |

- The site-specific nucleotide frequencies matrix shows that the number of occurrences of each of the four nucleotides at each location.

|  | sum | 1 | 2 | 3 | 4 | 5 | 6 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | 276 | 12 | 4 | 7 | 9 | 6 | 11 |
| C | 272 | 6 | 7 | 5 | 6 | 9 | 11 |
| G | 236 | 5 | 7 | 6 | 5 | 5 | 3 |
| T | 251 | 6 | 11 | 11 | 9 | 9 | 4 |

#### 2. Predictive Update.

- Firstly, we need **_N_** random numbers ranging from 1 to **_N_** and use these numbers as an index to choose the sequences sequentially to update the site-specific distribution of nucleotides and the associated frequencies. (Using a random series of numbers instead of choosing sequences according to the input order can decrease the likelihood of trapping Gibbs sampler within a local optimum.)
- Then we scan the sequence in a sliding window of length **_L_** to get position weight probability score (PPMS) for each window. 

![image.png](https://cdn.nlark.com/yuque/0/2021/png/2965461/1633510602157-a2d864c5-32ac-4c0e-8caa-0dd8f6a2885b.png#clientId=u9bd19277-e55c-4&from=paste&height=196&id=u3fca4f2b&margin=%5Bobject%20Object%5D&name=image.png&originHeight=392&originWidth=1174&originalType=binary&ratio=1&size=57452&status=done&style=none&taskId=ud64bae09-e319-4342-b58c-cc4e4c2ad74&width=587)​

- The following formula can be used to calculate PPMS. For each scanned sequence, we expect to get a new motif from this sequence that has the highest PPMS. The old ones in the motif set will be replaced.



![](https://cdn.nlark.com/yuque/__latex/422f46c3fd010a2263849613146f6307.svg#card=math&code=f_%7Bi%C2%B7pseudo%7D%3D%5Calpha%20f_i%5C%5C%0A%5C%5C%0Af_%7Bpseudo%7D%3D%5Csum%20f_%7Bi%C2%B7pseudo%7D%5C%5C%0A%5C%5C%0Ap_%7Bij%7D%3D%5Cfrac%20%7Bf_%7Bij%7D%2Bf_%7Bi%C2%B7pseudo%7D%7D%20%0A%7BN%2Bf_%7Bpseudo%7D%7D%5C%5C%0A%5C%5C%0A%5Ctext%20%7BPPMS%7D%3D%5Cprod_%7Bj%3D1%7D%5EL%20%5Cfrac%20%7Bp_%7Bij%7D%7D%20%7Bp_i%7D%20%5Ctext%7B%20%20%20%20%20%7D%20%7B_%7B%28i%5Ctext%7B%20%28A%20C%20G%20T%29%2C%20%20%20%7D%20j%5Ctext%7B%20%281%202%20...%20L%29%29%7D%7D%7D%5C%5C%20&id=hS0iD)

- We repeat this process for the rest of the sequences. After the last sequence has been updated, we have obtained a new set of motifs, together with the PPM. At this point we compute a weighted alignment score (**𝐹**) as follows: 

![](https://cdn.nlark.com/yuque/__latex/a84028e2ce1d9fdbd3d7fe71c4816c43.svg#card=math&code=F%3D%5Csum_%7Bi%3D1%7D%5EM%5Csum_%7Bj%3D1%7D%5ELC_%7Bij%7D%5Ctext%7BPPM%7D_%7Bij%7D%5C%5C&id=IDegs)

   - Each time when we get a new set of motifs, we compute a new **𝐹** value. **𝐹** is a measure of the quailty of alignment of the motifs. The larger the **𝐹** value, the better.
   - **_L_** is the motif width, and **_M_** is the number of different symbols in the sequences (4 for nucleotide and 20 for amino acid sequences).
- The predictive updating is repeated multiple times and previously stored locally optimal solutions are replaced by better ones. After going through all the sequences **1000** times, we finally find one of the locally optimal solutions.

![image.png](https://cdn.nlark.com/yuque/0/2021/png/2965461/1633528970680-cd182239-30dc-49f1-ac23-438c9096805a.png#clientId=u9bd19277-e55c-4&from=paste&height=426&id=uf7e71dc8&margin=%5Bobject%20Object%5D&name=image.png&originHeight=851&originWidth=2320&originalType=binary&ratio=1&size=111693&status=done&style=none&taskId=uec8b28da-a57a-41e1-8664-5a1b85654c3&width=1160)
#### 3. Convergence
Convergence is typically declared when two local solutions are identical.
#### 4. Output
Once Gibbs sampler has been finished, the output will be presented in six parts, including the following:

- Start and end time of the task.
- The final F value.
- The number of iterations.
- The motifs found where they start in the intron sequence.
- The site-specific nucleotide frequencies for the motifs found.
- The position weight matrix (PWM). 
```shell
Start: 2021-10-06 21:58
Gibbs sampler reached convergence after 10 iterations. 

Global alignment score (F) =  2799.25602709271 

Motif
     0  1  2  3  4  5  6
0   T  C  A  T  G  A  C
1   T  C  A  T  G  A  C
2   T  C  A  T  G  A  C
3   T  C  A  T  G  A  C
4   T  C  A  T  G  A  C
5   T  C  A  T  G  A  C
6   T  C  A  T  G  A  C
7   T  C  A  T  G  A  C
8   T  C  A  T  G  A  C
9   T  C  A  T  G  A  C
10  T  C  A  T  G  A  C
11  T  C  A  T  G  A  C
12  T  C  A  T  G  A  C
13  T  C  A  T  G  A  C
14  T  C  A  T  G  A  C
15  T  C  A  T  G  A  C
16  T  C  A  T  G  A  C
17  T  C  A  T  G  A  C
18  T  C  A  T  G  A  C
19  T  C  A  T  G  A  C
20  T  C  A  T  G  A  C
21  T  C  A  T  G  A  C
22  T  C  A  T  G  A  C
23  T  C  A  T  G  A  C
24  T  C  A  T  G  A  C
25  T  C  A  T  G  A  C
26  T  C  A  T  G  A  C
27  T  C  A  T  G  A  C
28  T  C  A  T  G  A  C
29  T  C  A  T  G  A  C
30  T  C  A  T  G  A  C
31  T  C  A  T  G  A  C
32  T  C  A  T  G  A  C
33  T  C  A  T  G  A  C
34  T  C  A  T  G  A  C
35  T  C  A  T  G  A  C
36  T  C  A  T  G  A  C
37  T  C  A  T  G  A  C
38  T  C  A  T  G  A  C
39  T  C  A  T  G  A  C
40  T  C  A  T  G  A  C
41  T  C  A  T  G  A  C
42  T  C  A  T  G  A  C
43  T  C  A  T  G  A  C
44  T  C  A  T  G  A  C
45  T  C  A  T  G  A  C
46  T  C  A  T  G  A  C
47  T  C  A  T  G  A  C
48  T  C  A  T  G  A  C
49  T  C  A  T  G  A  C
50  T  C  A  T  G  A  C
51  T  C  A  T  G  A  C
52  T  C  A  T  G  A  C
53  T  C  A  T  G  A  C
54  T  C  A  T  G  A  C
55  T  C  A  T  G  A  C
56  T  C  A  T  G  A  C
57  T  C  A  T  G  A  C
58  T  C  A  T  G  A  C
59  T  C  A  T  G  A  C
60  T  C  A  T  G  A  C
61  T  C  A  T  G  A  C
62  T  C  A  T  G  A  C
63  T  C  A  T  G  A  C
64  T  C  A  T  G  A  C
65  T  C  A  T  G  A  C
66  T  C  A  T  G  A  C
67  T  C  A  T  G  A  C
68  T  C  A  T  G  A  C
69  T  C  A  T  G  A  C
70  T  C  A  T  G  A  C
71  T  C  A  T  G  A  C
72  T  C  A  T  G  A  C
73  T  C  A  T  G  A  C
74  T  C  A  T  G  A  C
75  T  C  A  T  G  A  C
76  T  C  A  T  G  A  C
77  T  C  A  T  G  A  C
78  T  C  A  T  G  A  C
79  T  C  A  T  G  A  C
80  T  C  A  T  G  A  C
81  T  C  A  T  G  A  C
82  T  C  A  T  G  A  C
83  T  C  A  T  G  A  C
84  T  C  A  T  G  A  C
85  T  C  A  T  G  A  C
86  T  C  A  T  G  A  C
87  T  C  A  T  G  A  C
88  T  C  A  T  G  A  C
89  T  C  A  T  G  A  C
90  T  C  A  T  G  A  C
91  T  C  A  T  G  A  C
92  T  C  A  T  G  A  C
93  T  C  A  T  G  A  C
94  T  C  A  T  G  A  C
95  T  C  A  T  G  A  C
96  T  C  A  T  G  A  C
97  T  C  A  T  G  A  C
98  T  C  A  T  G  A  C
99  T  C  A  T  G  A  C 

Final site-specific counts:
      A    C    G    T
1    0    0    0  100
2    0  100    0    0
3  100    0    0    0
4    0    0    0  100
5    0    0  100    0
6  100    0    0    0
7    0  100    0    0 

Final PWM:
            A          C          G          T
1 -31.219281 -31.219281 -31.219281   2.000000
2 -31.219281   2.000000 -31.219281 -31.219281
3   2.000000 -31.219281 -31.219281 -31.219281
4 -31.219281 -31.219281 -31.219281   2.000000
5 -31.219281 -31.219281   2.000000 -31.219281
6   2.000000 -31.219281 -31.219281 -31.219281
7 -31.219281   2.000000 -31.219281 -31.219281 

Finish: 2021-10-06 21:59
```
It is possible that the input sequences may contain two or more different biologically significant motifs. This scripts will always end up with the strongest motif and miss all other motifs.
### Usage
```shell
usage: python gibbs.py -i read.fasta -L (int) [-o output/]

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        The path of fasta file, the suffix must be '.fa' or '.fasta'. [Required]
  -o OUTPUT, --output OUTPUT
                        output folder[Required]
  -L [LENGTH], --length [LENGTH]
                        The length of motif expected to get
```
