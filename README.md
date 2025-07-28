# 🏓 NICEgame exercise 🏓

This repo contains code and files used during a 2 week visit to the [EPFL LCSB](https://www.epfl.ch/labs/lcsb/). The goal was to get hands on experience using the [NICEgame](https://www.pnas.org/doi/10.1073/pnas.2211197119) workflow to curate metabolic models based on experimental data. Through these exercises we also learn about thermodynamic-based flux analysis (TFA), BridgIT, and ATLASx.

#### 🌡️ [Exercise 1](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_1): matTFA tutorial
#### 🐤 [Exercise 2](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_2): NICEgame tutorial
#### 💻 [Exercise 3](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_3): Essentiality prediction
#### ⚖️ [Exercise 4](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_4): Essentiality evaluation
#### 📑 [Exercise 5](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_5): Reproduction of result from NICEgame paper (AMAOTr)
#### 🧫 [Exercise 6](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_6): NICEgame with BIOLOG data
#### 🌉 [Exercise 7](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_7): Submitting a reaction query to the BridgIT server
#### 🗺️ [Exercise 8](https://github.com/franciscozorrilla/NICEgame_exercise/tree/master/Exercise_8): Submitting a metabolite query to the ATLASx server

## Requirements

- MATLAB: v2020b recommended for apple silicon
- CPLEX: v12.9 (full version free for academic use)

## Setup

After you have installed MATLAB and CPLEX, make sure that you add CPLEX to the MATLAB search path (accessed via the `Set path` button). Next, follow the installation instruction from the respective repos to install the following tools:

- [NICEgame](https://github.com/EPFL-LCSB/NICEgame)
- [matTFA](https://github.com/EPFL-LCSB/matTFA)

After cloning and pulling the repos, make sure to also add them in the MATLAB search path.

## License
The software in this repository is put under an APACHE licensing scheme

## Reference
If you use NICEgame (or any code in this repo) for your research, please consider citing the original publciation:

> E. Vayena, A. Chiappino-Pepe, H. MohammadiPeyhani, Y. Francioli, N. Hadadi, M. Ataman, J. Hafner, S. Pavlou, & V. Hatzimanikatis, A workflow for annotating the knowledge gaps in metabolic reconstructions using known and hypothetical reactions, Proc. Natl. Acad. Sci. U.S.A. 119 (46) e2211197119, https://doi.org/10.1073/pnas.2211197119 (2022).
