# Sparrow Search Algorithm

The project for the class `Swarm Intellgence.`

## Reference

> Xue, Jiankai, and Bo Shen. "A novel swarm intelligence optimization approach: sparrow search algorithm." Systems science & control engineering 8.1 (2020): 22-34.

>  [SSA python on github](https://github.com/changliang5811/SSA_python)


## Execution

```
g++ ./SSA.cpp -o ssa
./ssa [test_function_number] [version] [dimension] [run]
```


## Test Function number

| Number    | Function name |
| :-------: |:-------------:|
| 1         | Ackley        |
| 2         | Griewank      |
| 3         | Katsuuras     |
| 4         | Michalewicz   |
| 5         | Rosenbrock    |
| 6         | Schwefel2.26" |
| 7         | bent cigar    |
| 8         | zakharov      |
| 9         | happy cat     |
| 10        | rastrigin     |
| 11        | hgbat         |

If you input "0", the program will test all function without Katsuuras.


## Version

There are 7 mechanisms in this program and the version "0" is the original method.

| Number    | Mechanisms                                |
| :-------: |:-----------------------------------------:|
| 1         | time-variate st                           |
| 2         | time-variate sd                           |
| 3         | scrouger movement adapt levy flight       |
| 4         | aware danger movement adapt levy flight   |
| 5         | version 3 + version 4                     |
| 6         | version 1 + version 3                     |
| 7         | version 1 + version 4                     |


