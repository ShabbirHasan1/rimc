import matplotlib.pyplot as plt
from argparse import ArgumentParser
from itertools import product
from numpy.random import rand, randint
from numpy import ndarray, exp, round, zeros
import time


class Ensemble:
    arr: ndarray | None = None
    def __init__(self, dim: int, params: dict) -> None:
        self.dim = dim
        self.arr = zeros((self.dim, self.dim))
        self.params = params
    
    def init_ensemble(self) -> None:
        if self.arr is None:
            raise TypeError("self.arr is None and not an numpy.ndarray")
        arr_flat = self.arr.flatten()
        options = [-1, 1]
        for i, _ in enumerate(arr_flat):
            arr_flat[i] = options[int(round(rand()))]

        self.arr = arr_flat.reshape(self.dim, self.dim)

    def _calc_site_energy(self, i: int, j: int) -> float:
        coupling_const = self.params['coupling_const']

        s_k = self.arr[i][j]
        site_energy = 0.0

        nrows, ncols = self.arr.shape

        offsets = [-1, 1]

        for n in offsets:
            x = (i + n) % nrows
            y = (j + n) % nrows

            site_energy += self.arr[x][j]
            site_energy += self.arr[i][y]
        site_energy *= -coupling_const * s_k

        return site_energy

    def _calc_site_energy_diff(self, i: int, j: int) -> None:
        site_energy_diff = self._calc_site_energy(i, j)
        site_energy_diff *= -2.0

        return site_energy_diff

    def calc_state_energy(self) -> float:
        nrows, ncols = self.arr.shape

        state_energy = 0.0

        for i, j in product(range(nrows), range(ncols)):
            state_energy += self._calc_site_energy(i, j)
        return state_energy

    def step(self) -> None:
        beta = self.params['beta']

        nrows, ncols = self.arr.shape

        i = randint(nrows)
        j = randint(ncols)
        r = rand()

        accept_ratio = 0
        exponent = self._calc_site_energy_diff(i, j)

        if exponent <= 0:
            accept_ratio = 1
        else:
            exponent *= -beta
            accept_ratio = exp(exponent)

        if r < accept_ratio:
            self._update_site(i, j)
        
    def _update_site(self, i: int, j: int) -> None:
        spin = self.arr[i][j]

        match spin:
            case -1:
                self.arr[i][j] = +1
            case 1:
                self.arr[i][j] = -1
            case _:
                raise ValueError("Invalid value, cannot update state")


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-d", "--dimension", type=int, default=10,
                        help="Dimensions of the NxN ensemble")

    args = parser.parse_args()
    dim = args.dimension
    
    my_params = {
        'coupling_const': 1.0,
        'beta': 1.0,
        'mag_field': 1.0,
    }

    my_ensemble = Ensemble(dim, my_params)
    print(my_ensemble.arr)

    plt.ion()
    my_ensemble.init_ensemble()
    print(my_ensemble.arr)

    plt.figure(1)
    plt.imshow(my_ensemble.arr)

    t_start = time.time()
    for _ in range(my_ensemble.arr.size):
        my_ensemble.step()
    t_end = time.time()

    print(my_ensemble.arr)
    print(f"Total Duration: {t_end - t_start} s for {my_ensemble.arr.size} steps")

    plt.figure(2)
    plt.imshow(my_ensemble.arr)
