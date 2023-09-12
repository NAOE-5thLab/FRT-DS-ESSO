# FRT-DS-ESSO

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC_BY--NC_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

`FRT-DS-ESSO` is a repository that provides datasets measured in free running tests of Esso Osaka, intended for the purpose of data-driven ship maneuvering modeling. The datasets have been utilized for modeling or system identification in our research papers, specifically in the works authored by [Wakita (2022)](https://doi.org/10.1007/s00773-021-00867-1) and [Miyauchi (2022)](https://doi.org/10.1007/s00773-022-00889-3).



## Experiment Overview

### Subject ship (Scaled model ship of Esso Osaka)

The free-running model tests were carried out in the experimental pond known as Inukai pond, located at Osaka University. These tests utilized a model ship of [VLCC M.V. Esso Osaka](http://www.aukevisser.nl/inter-2/id427.htm). 

**An image taken at inukai pond**

![Photo of the subject model ship: M.V. ESSO OSAKA](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs00773-021-00867-1/MediaObjects/773_2021_867_Fig2_HTML.jpg?as=webp)

**Principal particular** 

| Item                                                         | Value  |
| :----------------------------------------------------------- | :----- |
| Length between perpendicular: $L_{\mathrm{pp}} \ (\mathrm{m})$ | 3.0    |
| Ship breadth: $B \ (\mathrm{m})$                             | 0.489  |
| Ship draft: $d \ (\mathrm{m})$                               | 0.201  |
| Diameter of propeller: $D_{\mathrm{p}} \ (\mathrm{m})$       | 0.084  |
| Area of Rudder: $A_{\mathrm{R}} \ \left(\mathrm{m}^2\right)$ | 0.0106 |
| Diameter of bow thruster: $D_{\mathrm{BT}} \ (\mathrm{m})$   | 0.050  |
| Diameter of stern thruster: $D_{\mathrm{ST}} \ (\mathrm{m})$ | 0.050  |
| Mass: $m \ (\mathrm{kg})$                                    | 244.6  |
| Longitudinal center of gravity: $x_{\mathrm{G}} \ (\mathrm{m})$ | 0.094  |
| Transverse projected area: $A_{\mathrm{T}} \ \left(\mathrm{m}^2\right)$ | 0.135  |
| Lateral projected area: $A_{\mathrm{L}} \ \left(\mathrm{m}^2\right)$ | 0.520  |
| Block coefficient: $C_{\mathrm{b}}$                          | 0.830  |

### Measurement Items

We employ two coordinate systems: the ship-fixed coordinate system denoted as $\mathrm{O}-xy$ and the space-fixed coordinate system denoted as $\mathrm{O}_0-x_0y_0$. These coordinate systems are defined as illustrated in the following right figure. Note that the origin of the $\mathrm{O}_0-x_0y_0$ system is located in close proximity to the corners of the berth.

<img src="./img/inukai.pdf" style="zoom:150%;" />

The observed variables have been outlined in the table provided below. These variables have been derived from measured data. For example, the ship position was observed by converting the GNSS receiver position to midship position. **If you require further information, please see [our paper](https://doi.org/10.1007/s00773-022-00889-3) for details on the measurement systems of the model ship.**

| Variable              | CSV header                     | Description                                 |
| --------------------- | ------------------------------ | ------------------------------------------- |
| $t$                   | `t [s]`                        | time of measurement.                        |
| $x_{0}$               | `x_position_mid [m]`           | midship position in $x_{0}$-axis direction. |
| $u$                   | `u_velo [m/s]`                 | longitudunal ship speed of midship.         |
| $y_{0}$               | `y_position_mid [m]`           | midship position in $y_{0}$-axis direction. |
| $v_{\mathrm{m}}$      | `vm_velo [m/s]`                | lateral ship speed of midship.              |
| $\psi$                | `psi_hat [rad]`                | ship heading angle.                         |
| $r$                   | `r_angvelo [rad/s]`            | yaw (heading) angular velocity.             |
| $n$                   | `n_prop [rps]`                 | revolution of propeller                     |
| $\delta$              | `delta_rudder [rad]`           | angle of rudder.                            |
| $U_{\mathrm{A}}$      | `wind_velo_relative_mid [m/s]` | relative wind velocity on midship.          |
| $\gamma_{\mathrm{A}}$ | `wind_dir_relative_mid [rad]`  | relative wind direction of midship.         |
| $U_{\mathrm{T}}$      | `wind_velo_true [m/s]`         | true wind velocity.                         |
| $\gamma_{\mathrm{T}}$ | `wind_dir_true [rad]`          | true wind direction.                        |

### Contents of free-running tests

The dataset contains measurement data from four types of free-running tests: turning test, zigzag test, random test, and berthing test.

#### Turning test

| File path                                        | $n \ \mathrm{[rps]}$ | $\delta \ \mathrm{[degree]}$ | File path                                        | $n \ \mathrm{[rps]}$ | $\delta \ \mathrm{[degree]}$ |
| ------------------------------------------------ | -------------------- | ---------------------------- | ------------------------------------------------ | -------------------- | ---------------------------- |
| `data/turning/turn_14-Oct-2020_14_08_51.csv`     | 8                    | 35                           | `data/turning/turn_14-Oct-2020_14_17_39.csv`     | 8                    | 35                           |
| `data/turning/turn_14-Oct-2020_14_29_49.csv`     | 10                   | 20                           | `data/turning/turn_14-Oct-2020_14_39_54.csv`     | 10                   | 20                           |
| `data/turning/turn_14-Oct-2020_14_56_07.csv`     | 10                   | -20                          | `data/turning/turn_14-Oct-2020_15_06_53.csv`     | 10                   | -20                          |
| `data/turning/turn_14-Oct-2020_15_21_57.csv`     | 10                   | -35                          | `data/turning/turn_14-Sep-2020_13_39_32.csv`     | 10                   | 35                           |
| `data/turning/turn_14-Sep-2020_13_51_45.csv`     | 10                   | 35                           | `data/turning/turn_14-Sep-2020_14_16_04.csv`     | 10                   | -35                          |
| `data/turning/turn_14-Sep-2020_14_50_41.csv`     | 8                    | 35                           | `data/turning/turn_14-Sep-2020_15_47_33.csv`     | 8                    | -35                          |
| `data/turning/turn_14-Sep-2020_15_58_08.csv`     | 10                   | 20                           | `data/turning/turn_cut_14-Sep-2020_15_58_08.csv` | 10                   | 20                           |
| `data/turning/turn_cut_14-Sep-2020_16_09_02.csv` | 10                   | -20                          |

#### Zigzag test

| File path                                     | $n \ \mathrm{[rps]}$ | $\delta \ \mathrm{[degree]}$ | File path                                     | $n \ \mathrm{[rps]}$ | $\delta \ \mathrm{[degree]}$ |
| --------------------------------------------- | -------------------- | ---------------------------- | --------------------------------------------- | -------------------- | ---------------------------- |
| `data/zigzag/zigzag_31-Jul-2020_13_04_24.csv` | 16.67                | $\pm$ 20                     | `data/zigzag/zigzag_31-Jul-2020_13_14_21.csv` | 16.67                | $\pm$ 30                     |
| `data/zigzag/zigzag_31-Jul-2020_13_22_52.csv` | 10                   | $\pm$ 15                     | `data/zigzag/zigzag_31-Jul-2020_13_29_19.csv` | 12                   | $\pm$ 15                     |
| `data/zigzag/zigzag_31-Jul-2020_13_42_53.csv` | 10                   | $\pm$ 30                     | `data/zigzag/zigzag_31-Jul-2020_13_50_28.csv` | 12                   | $\pm$ 30                     |
| `data/zigzag/zigzag_31-Jul-2020_13_57_45.csv` | 15                   | $\pm$ 20                     | `data/zigzag/zigzag_31-Jul-2020_14_03_39.csv` | 12                   | $\pm$ 20                     |
| `data/zigzag/zigzag_31-Jul-2020_14_10_05.csv` | 12                   | $\pm$ 20                     |


#### Random test

The random test aims to contain all possible values of control inputs and states to reduce necessary data for training and tests. In this test, control inputs were selected randomly by the human operator to collect datasets including various ship and actuator states. Please see [our paper](https://doi.org/10.1007/s00773-022-00889-3) for details of this tests.

#### Berthing test

The berthing maneuver was conducted in the center of the pond without the berth wall and controlled manually by the operator.



## Correspondence table with the paper

- [Wakita (2022)](https://doi.org/10.1007/s00773-021-00867-1)

| File path                                         | Train-TZB | Train-TZRB | Train-TZRB+ | Test (TZRB) |
| ------------------------------------------------- | :-------: | :--------: | :---------: | :---------: |
| `data/random/random_03-Aug-2020_11_04_37.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/random/random_06-Aug-2020_16_24_25_1.csv`   |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/random/random_06-Aug-2020_16_24_25_3.csv`   |     ✖️     |     ✖️      |      ✖️      |      ⭕️      |
| `data/random/random_06-Aug-2020_16_24_25_4.csv`   |     ✖️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/random/random_06-Aug-2020_16_24_25_5.csv`   |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/random/random_06-Aug-2020_16_24_25_6.csv`   |     ✖️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/random/random_31-Jul-2020_12_29_58.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/random/random_31-Jul-2020_12_48_39.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/random/random_31-Jul-2020_14_52_18.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Oct-2020_14_08_51.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Oct-2020_14_17_39.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Oct-2020_14_29_49.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Oct-2020_14_39_54.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Oct-2020_14_56_07.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Oct-2020_15_06_53.csv`      |     ✖️     |     ✖️      |      ✖️      |      ⭕️      |
| `data/turning/turn_14-Oct-2020_15_21_57.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Sep-2020_13_39_32.csv`      |     ⭕️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Sep-2020_13_51_45.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Sep-2020_14_16_04.csv`      |     ⭕️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Sep-2020_14_50_41.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Sep-2020_15_47_33.csv`      |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_14-Sep-2020_15_58_08.csv`      |     ⭕️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/turning/turn_cut_14-Sep-2020_15_58_08.csv`  |     ✖️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/turning/turn_cut_14-Sep-2020_16_09_02.csv`  |     ⭕️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_13_04_24.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_13_14_21.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_13_22_52.csv`     |     ⭕️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_13_29_19.csv`     |     ⭕️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_13_42_53.csv`     |     ✖️     |     ✖️      |      ✖️      |      ⭕️      |
| `data/zigzag/zigzag_31-Jul-2020_13_50_28.csv`     |     ⭕️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_13_57_45.csv`     |     ⭕️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_14_03_39.csv`     |     ⭕️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/zigzag/zigzag_31-Jul-2020_14_10_05.csv`     |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/berthing/berthing_12-Nov-2020_11_07_23.csv` |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/berthing/berthing_12-Nov-2020_11_14_47.csv` |     ⭕️     |     ⭕️      |      ⭕️      |      ✖️      |
| `data/berthing/berthing_12-Nov-2020_11_21_46.csv` |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/berthing/berthing_12-Nov-2020_11_38_50.csv` |     ✖️     |     ✖️      |      ✖️      |      ⭕️      |
| `data/berthing/berthing_12-Nov-2020_11_49_21.csv` |     ✖️     |     ✖️      |      ✖️      |      ⭕️      |
| `data/berthing/berthing_12-Nov-2020_13_09_11.csv` |     ✖️     |     ✖️      |      ⭕️      |      ✖️      |
| `data/berthing/berthing_12-Nov-2020_13_21_03.csv` |     ✖️     |     ✖️      |      ✖️      |      ⭕️      |
| `data/berthing/berthing_12-Nov-2020_13_37_05.csv` |     ⭕️     |     ⭕️      |      ⭕️      |      ✖️      |

- [Miyauchi (2022)](https://doi.org/10.1007/s00773-022-00889-3)

| File path                                             | Train-R | Train-TR | Train-TRZ | Test |
| ----------------------------------------------------- | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| `data/random/random_03-Aug-2020_11_04_37.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/random/random_06-Aug-2020_16_24_25_1.csv`   |                              ⭕️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/random/random_06-Aug-2020_16_24_25_3.csv`   |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ⭕️                               |
| `data/random/random_06-Aug-2020_16_24_25_4.csv`   |                              ⭕️                               |                              ⭕️                               |                              ⭕️                               |                              ✖️                               |
| `data/random/random_06-Aug-2020_16_24_25_5.csv`   |                              ⭕️                               |                              ⭕️                               |                              ⭕️                               |                              ✖️                               |
| `data/random/random_06-Aug-2020_16_24_25_6.csv`   |                              ⭕️                               |                              ⭕️                               |                              ⭕️                               |                              ✖️                               |
| `data/random/random_31-Jul-2020_12_29_58.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/random/random_31-Jul-2020_12_48_39.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/random/random_31-Jul-2020_14_52_18.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_14_08_51.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_14_17_39.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_14_29_49.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_14_39_54.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_14_56_07.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_15_06_53.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Oct-2020_15_21_57.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Sep-2020_13_39_32.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Sep-2020_13_51_45.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Sep-2020_14_16_04.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Sep-2020_14_50_41.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ⭕️                               |
| `data/turning/turn_14-Sep-2020_15_47_33.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_14-Sep-2020_15_58_08.csv`      |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_cut_14-Sep-2020_15_58_08.csv`  |                              ✖️                               |                              ⭕️                               |                              ✖️                               |                              ✖️                               |
| `data/turning/turn_cut_14-Sep-2020_16_09_02.csv`  |                              ✖️                               |                              ⭕️                               |                              ⭕️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_04_24.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_14_21.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_22_52.csv`     |                              ✖️                               |                              ✖️                               |                              ⭕️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_29_19.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_42_53.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_50_28.csv`     |                              ✖️                               |                              ✖️                               |                              ⭕️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_13_57_45.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_14_03_39.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/zigzag/zigzag_31-Jul-2020_14_10_05.csv`     |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ⭕️                               |
| `data/berthing/berthing_12-Nov-2020_11_07_23.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/berthing/berthing_12-Nov-2020_11_14_47.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/berthing/berthing_12-Nov-2020_11_21_46.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/berthing/berthing_12-Nov-2020_11_38_50.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/berthing/berthing_12-Nov-2020_11_49_21.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ⭕️                               |
| `data/berthing/berthing_12-Nov-2020_13_09_11.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |
| `data/berthing/berthing_12-Nov-2020_13_21_03.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ⭕️                               |
| `data/berthing/berthing_12-Nov-2020_13_37_05.csv` |                              ✖️                               |                              ✖️                               |                              ✖️                               |                              ✖️                               |



## Citation

If you use dataset in a scientific publication, we would appreciate references to the following papers:

- **Wakita, K., Maki, A., Umeda, N., Miyauchi, Y., Shimoji, T., Rachman, D. M., & Akimoto, Y. (2022). On neural network identification for low-speed ship maneuvering model. *Journal of Marine Science and Technology*, *27*(1), 772–785. https://doi.org/10.1007/s00773-021-00867-1**
- **Miyauchi, Y., Maki, A., Umeda, N., Rachman, D. M., & Akimoto, Y. (2022). System parameter exploration of ship maneuvering model for automatic docking/berthing using CMA-ES. *Journal of Marine Science and Technology*, *27*(2), 1065–1083. https://doi.org/10.1007/s00773-022-00889-3**

Biblatex entry:

```
@article{Wakita2022,
   author = {Kouki Wakita and Atsuo Maki and Naoya Umeda and Yoshiki Miyauchi and Tohga Shimoji and Dimas M Rachman and Youhei Akimoto},
   doi = {10.1007/s00773-021-00867-1},
   issn = {1437-8213},
   issue = {1},
   journal = {Journal of Marine Science and Technology},
   pages = {772-785},
   title = {On neural network identification for low-speed ship maneuvering model},
   volume = {27},
   url = {https://doi.org/10.1007/s00773-021-00867-1},
   year = {2022},
}
```

```
@article{Miyauchi2022,
   author = {Yoshiki Miyauchi and Atsuo Maki and Naoya Umeda and Dimas M Rachman and Youhei Akimoto},
   doi = {10.1007/s00773-022-00889-3},
   issn = {1437-8213},
   issue = {2},
   journal = {Journal of Marine Science and Technology},
   pages = {1065-1083},
   title = {System parameter exploration of ship maneuvering model for automatic docking/berthing using CMA-ES},
   volume = {27},
   url = {https://doi.org/10.1007/s00773-022-00889-3},
   year = {2022},
}
```



## License

The dataset is released under the [Creative Commons BY-NC 4.0 License](https://creativecommons.org/licenses/by-nc/4.0/).

[![CC BY-NC 4.0](https://mirrors.creativecommons.org/presskit/buttons/88x31/png/by-nc.png)](https://creativecommons.org/licenses/by-nc/4.0/)



