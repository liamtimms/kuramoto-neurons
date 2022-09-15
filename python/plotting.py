import os

import imageio
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def load_data(filename):
    return np.load(filename)


def plot_phi_vs_t(save_name, data):
    """TODO: Docstring for plot_all.

    :arg1: TODO
    :returns: TODO

    """
    flat_data = data.reshape(-1, data.shape[-1])
    plt.plot(flat_data.T)
    plt.savefig(save_name)
    plt.close()
    return


def plot_1D(input_dir, identifier):
    """TODO: Docstring for plot_1D.

    :arg1: TODO
    :returns: TODO

    """
    phi_load_name = os.path.join(input_dir, f"phi_{identifier}.npy")
    data = load_data(phi_load_name)
    save_name = os.path.join(input_dir, f"phi-vs-t_{identifier}.png")
    plot_phi_vs_t(save_name, data)
    filenames = []
    for t in tqdm(range(data.shape[-1])):
        # don't plot every frame
        # step = data.shape[-1] / 20
        step = 20
        if (t % step) == 0:
            row = data[:, t].ravel() % (2 * np.pi)  # ravel to make it 1D
            frequency = (data[:, t].ravel() - data[:, t - step].ravel()) / step
            for phi, omega in zip(row, frequency):
                plt.polar(phi, omega, "ro")
                # create file name and append it to a list
            filename = f"{t}.png"
            filenames.append(filename)

            # save frame
            plt.savefig(filename)
            plt.close()  # build gif

    save_name = os.path.join(input_dir, f"phi_{identifier}.gif")
    with imageio.get_writer(save_name, mode="I") as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)

    pass


def plot_2D(input_dir, identifier):
    """TODO: Docstring for plot_2D.

    :arg1: TODO
    :returns: TODO

    """
    phi_load_name = os.path.join(input_dir, f"phi_{identifier}.npy")
    data = load_data(phi_load_name)
    save_name = os.path.join(input_dir, f"phi-vs-t_{identifier}.png")
    plot_phi_vs_t(save_name, data)
    filenames = []
    for t in tqdm(range(data.shape[-1])):
        # don't plot every frame
        # step = data.shape[-1] / 20
        step = 20
        if (t % step) == 0:
            slice = data[:, :, t].squeeze() % (2 * np.pi)  # squeeze to 2D
            plt.imshow(slice, vmin=0, vmax=2 * np.pi, cmap="hsv")
            # create file name and append it to a list
            filename = f"{t}.png"
            filenames.append(filename)

            # save frame
            plt.savefig(filename)
            plt.close()  # build gif

    save_name = os.path.join(input_dir, f"phi_{identifier}.gif")
    make_gif(save_name, filenames)

    pass


def make_gif(save_name, filenames):
    with imageio.get_writer(save_name, mode="I") as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)


def main():
    input_dir = "../kuramoto-rust/target/release/test_output"

    # parameters of run
    dimension = 1
    n = 20
    timesim = 6000
    timemetric = 2
    spreadinomega = 0.25
    g = 1
    epsilon = 0.1
    clustersize = 0
    drivingfrequency = 1

    identifier = (f"N{n}-D{dimension}-t{timesim}-T{timemetric}" +
                  f"-ostd{spreadinomega}-g{g}-e{epsilon}" +
                  f"-C{clustersize}-do{drivingfrequency}")

    if dimension == 1:
        plot_1D(input_dir, identifier)
    elif dimension == 2:
        plot_2D(input_dir, identifier)
    else:
        print(f"Dimension {dimension} not implemented")


if __name__ == "__main__":
    main()
