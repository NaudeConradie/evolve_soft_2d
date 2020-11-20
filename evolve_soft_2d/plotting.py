##  Plotting functions

#   Imports
import math
import matplotlib.pyplot as plot
import pandas
import seaborn

from evolve_soft_2d import file_paths

################################################################################

def histogram(
    template,
    tm: str,
    data: list,
    t: str,
    y: str,
    x: str,
    bins = "auto",
    color: str = "b",
    ) -> None:
    """Plot a histogram

    Parameters
    ----------
    template
        The unit template parameters
    tm : str
        The timestamp of the current simulation
    data : list
        The data to be plotted
    t : str
        The title of the graph
    y : str
        The label of the y-axis
    x : str
        The label of the x-axis
    bins : optional
        The bin settings, by default "auto"
    color : str, optional
        The colour of the graph, by default "b"
    """

    #   Determine the maximum x-axis value of the graph
    bin_max = math.ceil(max(data))

    #   Open a figure
    plot.figure()

    #   Plot the histogram
    plot.rcParams.update({"figure.figsize":(7, 5), "figure.dpi":100})
    plot.hist(data, bins = bins, color = color)
    plot.gca().set(title = t, ylabel = y, xlabel = x)
    plot.xlim(0, bin_max)

    # #   Show the plot
    # plot.show()

    #   Save the figure
    save_plot(template, t, tm)

    return

################################################################################

def hist_all(
    template,
    tm,
    data: pandas.DataFrame,
    ) -> None:

    data_col = [i for i in data.columns]

    for i in data_col:

        seaborn.displot(data, x = i)
        # seaborn.displot(data, x = i, bins = 20)

        #   Save the figure
        save_plot(template, i, tm)

    return

################################################################################

def hist(
    template,
    tm,
    data: pandas.DataFrame,
    x: str,
    ) -> None:

    seaborn.displot(data, x = x)

    #   Save the figure
    save_plot(template, x, tm)

    return

################################################################################

def boxp(
    template,
    tm,
    data: pandas.DataFrame,
    x: str,
    y: str,
    ) -> None:

    seaborn.boxplot(x = x, y = y, data = data)

    #   Save the figure
    save_plot(template, x + "_vs_" + y, tm)

    return

################################################################################

def boxp_all(
    template,
    tm,
    data: pandas.DataFrame,
    x: list,
    y: list,
    ) -> None:

    for i in x:

        for j in y:

            seaborn.boxplot(x = i, y = j, data = data)

            #   Save the figure
            save_plot(template, i + "_vs_" + j, tm)

    return

################################################################################

def boxp_melt(
    template,
    tm,
    data: pandas.DataFrame,
    x: list,
    y: list,
    ) -> None:

    for i in x:

        data_melt = data.melt(id_vars = i, value_vars = y, var_name = "Type of Energy", value_name = "Energy (mJ)") 

        seaborn.boxplot(x = i, y = "Energy (mJ)", hue = "Type of Energy", data = data_melt)

        #   Save the figure
        save_plot(template, i + "_vs_" + "Energy (mJ)", tm)

    return

################################################################################

def scat_all(
    template,
    tm,
    data: pandas.DataFrame,
    x: list,
    y: list,
    ) -> None:

    for i in x:

        for j in y:

            seaborn.relplot(x = i, y = j, data = data)

            #   Save the figure
            save_plot(template, i + "_vs_" + j, tm)

    return

################################################################################

def scatterplot(
    template,
    tm: str,
    x_data: list,
    y_data: list,
    t: str,
    x_l: str,
    y_l: str,
    color: str = "b",
    marker: str = "o"
    ) -> None:
    """Plot a scatter plot

    Parameters
    ----------
    template
        The unit template parameters
    tm : str
        The timestamp of the current simulation
    x_data : list
        The data to be plotted on the x-axis
    y_data : list
        The data to be plotted on the y-axis
    t : str
        The title of the graph
    y : str
        The label of the y-axis
    x : str
        The label of the x-axis
    color : str, optional
        The colour of the graph, by default "b"
    marker : str, optional
        The plot markers, by default "o"
    """

    #   Determine the maximum x-axis value of the graph
    x_max = math.ceil(max(x_data))

    #   Open a figure
    plot.figure()

    #   Plot the scatter plot
    plot.rcParams.update({"figure.figsize":(7, 5), "figure.dpi":100})
    plot.scatter(x_data, y_data, c = color, marker = marker)
    plot.gca().set(title = t, ylabel = y_l, xlabel = x_l)
    plot.xlim(0, x_max)

    # #   Show the plot
    # plot.show()

    #   Save the figure
    save_plot(template, t, tm)

    return

################################################################################

def lreg_all(
    template,
    tm,
    data: pandas.DataFrame,
    x: list,
    y: list,
    order: int = 1
    ) -> None:

    for i in x:

        for j in y:

            seaborn.regplot(x = i, y = j, data = data, order = order)

            #   Save the figure
            save_plot(template, i + "_vs_" + j, tm)

    return

################################################################################

def grid_plot(
    template,
    tm,
    data: pandas.DataFrame,
    l: str,
    ) -> None:

    seaborn.heatmap(data, linewidths=1, linecolor='black', cmap = "Greys")

    #   Save the figure
    save_plot(template, l, tm)

    return

################################################################################

def plot_all(
    template,
    v: list,
    n_e: list,
    l: list,
    tm: str,
    ) -> None:
    """Plot all desired figures

    Parameters
    ----------
    template
        The unit template parameters
    v : list
        The data to be plotted
    n_e : list
        The list of the number of elements removed from every element
    l : list
        The list of labels of the data
    tm : str
        The timestamp of the current simulation
    """

    scatterplot(template, tm, v[0], v[1], "Constraint Energy X vs Y", "Constraint Energy X (J)", "Constraint Energy Y (J)")
    scatterplot(template, tm, v[3], v[4], "Internal Energy X vs Y", "Internal Energy X (J)", "Internal Energy Y (J)")

    scatterplot(template, tm, n_e, v[6], "Elements Removed vs Hausdorff Distance", "Number of Elements Removed", "Hausdorff Distance")

    # #   Loop through the types of data
    # for i in range(0, len(v)):

    #     #   Plot the histogram
    #     histogram(template, tm, v[i], l[i], "Frequency", "Energy (J)")

    #     #   Plot the scatterplot
    #     scatterplot(template, tm, n_e, v[i], l[i], "Energy (J)", "Number of Elements Removed")

    # #   Plot a scatterplot
    # scatterplot(template, tm, v[0], v[1], "Constraint Energy (J)", "Y-direction", "X-direction")

    return

################################################################################

def save_plot(
    template,
    t: str,
    tm: str,
    ) -> None:
    """Save a figure

    Parameters
    ----------
    template
        The unit template parameters
    t : str
        The title of the graph
    tm : str
        The timestamp of the current simulation
    """    

    #   Create the file path of the figure
    fp_p = file_paths.create_fp_file(template, t + tm, "g")

    plot.tight_layout()

    #   Save the figure
    plot.savefig(fp_p, dpi = 300)

    #   Close the figure
    plot.close()

    return