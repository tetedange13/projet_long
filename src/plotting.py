#!/usr/bin/env python3

"""
Module gethering all functions requiered to generate the different plots
"""


import matplotlib.pyplot as plt
import matplotlib.ticker as tck


def display_curve(levels_x, levels_x_rev, TMscores_y, TMscores_y_rev,
                  res_gdt, res_gdt_rev,
                  ref_pdb_id, peeled_pdb_id, parMATT_TMscore, ref_TMscore):
    """
    Display the curve of the values of TMscore (between aligned PU and reference
    pdb), according to the level considered

    Args:
        levels_x: Levels considered (list)
        TMscores_y: Values of TMscore (list)
        pdb_ref_name: Name (str) of the reference pdb
    """
    fig = plt.figure()
    axis = plt.subplot(111)
    plt.title("TMscore according to the level considered")

    # A list of plot to add different scores in the same plot
    plots = []
    plots.append(plt.plot(levels_x, TMscores_y, 'bx-'))
    plots.append(plt.plot(levels_x_rev, TMscores_y_rev, 'gx-'))

    # Get maximum number of PUs:
    max_nb_PU = max(max(levels_x), max(levels_x_rev))
    # Display horizontal line for parMATT TMscore:
    plots.append(plt.plot([0, max_nb_PU],
                          [parMATT_TMscore] * 2, 'k-'))
    # Display horizontal line for reference (normal) TMscore:
    plots.append(plt.plot([0, max_nb_PU],
                          [ref_TMscore] * 2, 'm-'))

    # gdt results:
    plots.append(plt.plot(levels_x, res_gdt, 'cx-'))
    plots.append(plt.plot(levels_x_rev, res_gdt_rev, 'rx-'))

    # Add legends and axis labels:
    plt.legend([plot[0] for plot in plots],
               (peeled_pdb_id + " peeled", ref_pdb_id + " peeled",
                "parMTT TMscore", 'TMscore ref'))
    plt.ylabel('TMscores between ' + ref_pdb_id + ' and ' + peeled_pdb_id)
    plt.xlabel('Number of PUs at each level of cutting')
    axis.xaxis.set_major_locator(tck.MaxNLocator(integer=True))

    fig.savefig("curves.pdf")
    plt.close(fig)


def disp_barplot(all_means):
    """
    Display a barplot with mean values of TMscores

    all_means: A dataFrame (pandas) with the values of mean for each score
    """
    fig = plt.figure()
    axis = plt.subplot(111)
    plt.title("Mean values of TMscore according to the method considered")


    barlist = axis.bar(all_means.index, all_means.values)
    fig.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=45)
    for bar in barlist:
        bar.set_color(list(np.random.random(size=3)))

    fig.savefig("barplot.pdf")
    plt.close(fig)
