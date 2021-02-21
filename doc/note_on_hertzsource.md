# Note on Hertz source implementation
2021.02.21 Kurama Okubo

This notebook describes the theory of source time function of ball drop test following

McLaskey, G.C., Glaser, S.D. Acoustic Emission Sensor Calibration for Absolute Source Measurements. J Nondestruct Eval 31, 157–168 (2012). https://doi.org/10.1007/s10921-012-0131-2

and its implementation in OpenSWPC.

## Theory of Hertzian source

Hertzian theory gives analytical solution of ball impact. From the equation 7 in McLaskey and Glaser (2012), it is written as follows:

$$
f(t) = \begin{cases}
f_{max} \left\{\sin (\pi t / t_c) \right\} ^ {3/2} & 0 \leq |t| \leq t_c,\\
0 & \text{otherwise,}\\
\end{cases}
$$

where the contact time $t_c$ is

$$
t_c = 4.53( 4\rho_1 \pi (\delta_1 + \delta_2) /3 )^{2/5} R_1 v_0 ^{-1/5}.
$$

$\rho_1$ is the density, $R_1$ is the radius and $v_0$ is the velocity of incoming ball.

The maximum force $ f_{max} $ is

$$
f_{max} = 1.917\rho_1 ^{3/5} (\delta_1 + \delta_2) ^{-2/5} R_1 ^2 v_0 ^{6/5},
$$

$\delta_i$ is defined as

$$
\delta_i = (1-\mu_i^2)(\pi E_i),
$$

where $E$ and $\mu$ are the Young's modulus and Poisson's ratio, respectively. The subscripts 1 and 2 indicate the material properties of ball and the base specimen, respectively.

## Typical values of constants

| property | description  | value |
|:------:|:------| :------|
|  $\rho_1$| density of ball (sus) | 7930 ^\*2 [kg/m^3] |
|  $R_1$  | radius of ball| 1.5e-3 [m] |
|  $v_0$  | incoming velocity from 50cm hight |3.13156 [m/s]|
|  $E_1$  | Young's modulus of ball| 197 ^\*3  [GPa]|
|  $\mu_1$| Poisson's ratio of ball| 0.3 ^\*3 |
|  $E_2$  |Young's modulus of rock specimen| 103 ^\*1 [GPa]|
|  $\mu_2$|Poisson's ratio of rock specimen | 0.31 ^\*1 |

*1 Fukuyama et al. (2016)

*2 [link](http://www.jssc.or.jp/ssba/structure/q&a/q&a.html#:~:text=%E4%B8%80%E8%88%AC%E9%8B%BC%E3%81%AE%E3%81%9B%E3%82%93%E6%96%AD%E5%BC%BE%E6%80%A7,%E3%81%BE%E3%81%9F%E3%80%81%E9%8B%BC%E7%A8%AE%E3%81%AB%E3%82%88%E3%82%8A%E7%95%B0%E3%81%AA%E3%82%8A%E3%81%BE%E3%81%99%E3%80%82&text=%E3%82%B9%E3%83%86%E3%83%B3%E3%83%AC%E3%82%B9%E9%8B%BC%E3%81%AE%E9%83%A8%E6%9D%90%E3%82%84,%E4%B8%80%E8%88%AC%E9%8B%BC%E3%81%A8%E7%95%B0%E3%81%AA%E3%82%8A%E3%81%BE%E3%81%99%E3%80%82)

*3 [link](https://d-engineer.com/cae/material.html)


## Strategy of implementation

In OpenSWPC, all source time functions (STF) are normalized so that $\int_0^\infty \dot{M}(t) dt = 1$. Computations are performed using this normalized STF for the sake of numerical stability. Then, when output the results, the magnitude is corrected with input M0.

To follow this manner, it is ideal that we normalize STF of Hertz source model by its total impulse. However, the integration is tricky as we need complicated integral form. While I just make a mathematical form of normalization factor for Hertz STF, I decided that we don't normalize by M0 (i.e. M0 = $f_z$ = 1.0 in code), which is enough to achieve the numerical stability based on the analysis of normalization factor described in the next section.

So the message is that **do not use Hertz STF for moment tensor source with M0 input** because the result does not reflect the true M0 because $\int_0^\infty f(t) dt \neq 1$.

Use Hertz STF for the body force mode with constants associated with dropping ball.

### Mathematical formulation of impulse for Hertz source
Our goal here is to compute

$$\int_{0}^{t_c} f(t) dt = \int_{0}^{t_c} f_{max} \left\{\sin (\pi t / t_c) \right\} ^ {3/2} dt,$$
which can be used as normalization factor equivalent to M0 in OpenSWPC.

$$\int_{0}^{t_c} f_{max} \left\{\sin (\pi t / t_c) \right\} ^ {3/2} dt = f_{max} \frac{t_c}{\pi} \int_{0}^{\pi} \left\{\sin (\phi) \right\} ^ {3/2} d\phi $$

Let
$$ I = \int_{0}^{\pi} \left\{\sin (\phi) \right\} ^ {3/2} d\phi, $$

Then
$$ \int_{0}^{t_c} f_{max} \left\{\sin (\pi t / t_c) \right\} ^ {3/2} dt = f_{max} \frac{t_c}{\pi}  I $$


$$I=\int_{0}^{\pi} \frac{\cos ^2 {\phi}} {2\sqrt{\sin{\phi}}} d\phi$$
$$=\int_{0}^{\pi} \frac{1}{2\sqrt{\sin{\phi}}} d\phi - \frac{1}{2}I$$

Thus
$$I=\frac{1}{3} \int_{0}^{\pi}\frac{1}{\sqrt{\sin{\phi}}}d\phi. $$

I couldn't find the value above, so computed numerically and found a value as follows:

$$\frac{1}{3} \int_{0}^{\pi}\frac{1}{\sqrt{\sin{\phi}}}d\phi=\frac{1}{3}*5.244115108583619$$

Thus
$$I = 1.748038369527873$$

Therefore, the normalization factor $k$ is written as:

$$ k = 1.748 f_{max} \dfrac{t_c}{\pi}.$$

Note that if using the typical values, $t_c \simeq 1e-5$ and $f_{max} = 131.0$. Then $k = 7.29e-4$. So the range of magnitude in the variables associated with source varies in the order of 1e-4, which is fine in terms of overflow of floating point when using double precision.

## Tips for code implementation

`m_source.F90` l218: $$M0 = \sqrt{f_x^2 +f_y^2 + f_z^2}$$ for body force mode. This is not used for ball drop as we normalize with $f_{max}$

`m_source.F90` l335: $f_z$ is normalized by 'M0', which is corrected when output in `m_output.F90`.


##　メモ
OpenSWPCでは桁落ち防止のためにSTFをM0で正規化して計算し、出力する段階でM0をかけて補正している。Hertzの震源関数を用いる場合でも正規化するのが理想的であるが、その積分値を解析的に求めるのが難しいため、また実装が複雑になるため本コードでは正規化せずに計算を行うことにした。その際、桁落ちに関して検討したところ、数値積分も利用して上記のように正規化の定数を推定したところおおよそ7.29e-4であった。つまりSTFを正規化する場合よりも5桁ほど絶対値が小さい状態で計算することになる。これは倍精度の桁数約16桁よりも十分小さいので正規化をしないことによる桁落ちの心配はないと結論づけた。
