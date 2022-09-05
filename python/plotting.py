import os

import imageio
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def load_data(filename):
    return np.load(filename)


def main():
    """TODO: Docstring for main.

    :arg1: TODO
    :returns: TODO

    """
    data = load_data(
        "../kuramoto-rust/target/release/phi_N20-time10000-D1.npy")
    filenames = []
    for t in tqdm(range(data.shape[-1])):
        # don't plot every frame
        if t % 20 == 0:
            row = data[:, t].ravel() % (2 * np.pi)
            frequency = (data[:, t].ravel() - data[:, t - 10].ravel()) / 10
            for phi, omega in zip(row, frequency):
                plt.polar(phi, omega, "ro")
                # create file name and append it to a list
            filename = f"{t}.png"
            filenames.append(filename)

            # save frame
            plt.savefig(filename)
            plt.close()  # build gif

    with imageio.get_writer("mygif.gif", mode="I") as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)

        plt.show()

    pass


if __name__ == "__main__":
    main()
