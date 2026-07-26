# Daily Weather data of Kumasi from 1960 to the present

A dataset containing Rain, maximum and minimum temperature, wind speed,
sunshine hours, and minimum and maximum Relative Humidity for Kumasi,
Ghana. The record from 1960 to 2015 is the historical station series
(see
[`kumasi1960_2015`](https://gowusu.github.io/sebkc/reference/kumasi1960_2015.md));
from 2016 onward it is extended with NASA POWER daily data retrieved
through [`weather`](https://gowusu.github.io/sebkc/reference/weather.md)
(latitude 6.72, longitude -1.60). For the appended period the incoming
shortwave radiation Rs is converted to sunshine hours `n` with the
inverse Angstrom-Prescott relation so that the whole series keeps a
single, homogeneous set of columns and can be passed to
[`ETo`](https://gowusu.github.io/sebkc/reference/ETo.md) and
[`kc`](https://gowusu.github.io/sebkc/reference/kc.md) without special
handling. Regenerate or extend it at any time by re-running
[`weather()`](https://gowusu.github.io/sebkc/reference/weather.md) for
the new dates and re-appending.

## Format

A data frame with 8 variables:

- DOY:

  Calendar date \[YYYY-mm-dd\]

- Tmax:

  Maximum Temperature \[oC\]

- Tmin:

  Minimum Temperature \[oC\]

- RHmin:

  Minimum Relative Humidity \[per cent\]

- RHmax:

  Maximum Relative Humidity \[per cent\]

- n:

  sunshine hours

- uz:

  wind speed \[m/s\]

- P:

  Daily precipitation \[mm\]

## Author

George Owusu
