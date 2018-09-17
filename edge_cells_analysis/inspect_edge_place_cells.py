import numpy as np
import os
import cPickle

from edge_vs_place_recurrence import find_edge_place_cells, plot_matrix
from bambi_playback_analysis.code.tools import plot_cell_activity_in_session
from cell_activity_in_different_environments import *
from zivlab.analysis import population
from bambi_playback_analysis.code import tools

PLACE_CELLS_PATH = r"D:\dev\edge_cells_analysis\place cells\two environments\1000 permutations"
OUTPUT_PATH = r"D:\dev\edge_cells_analysis\place cells\two environments\1000 permutations\velocity threshold for edge cells figures"

def load_place_cells(filename):
    place_cells_data = cPickle.load(open(filename))
    sessions_place_cells = place_cells_data['sessions_place_cells']
    env_a_place_cells = [x[0][0] for x in sessions_place_cells]
    env_b_place_cells = [x[1][0] for x in sessions_place_cells]
    velocity_threshold = place_cells_data['velocity_threshold']
    min_number_of_events = place_cells_data['min_number_of_events']

    return env_a_place_cells, env_b_place_cells, velocity_threshold, min_number_of_events


def main():
    for mouse in MICE:
        print mouse
        place_cells_filename = \
            os.path.join(PLACE_CELLS_PATH, "%s_place_cells.pkl" %mouse)
        env_a_place_cells, env_b_place_cells, velocity_threshold, min_number_of_events = \
            load_place_cells(place_cells_filename)
        all_days_overlap = []
        all_days_recurrence = []
        for day in DAYS:
            print "day %d" % day
            try:
                # Load session data
                day_path = os.path.join(DATA_PATH, mouse, 'day%d' % day)
                env_a_path = os.path.join(day_path, 'envA')
                env_b_path = os.path.join(day_path, 'envB')
                cell_reg_filename = os.path.join(
                    CELL_REGISTRATION_PATH, mouse, "cellRegisteredFixed.mat")
                cell_registration = matlab.load_cell_registration(
                    cell_reg_filename)
                events_a, movement_a = load_two_env_session(
                    env_a_path, session_ind=2*(day-1),
                    cell_registration=cell_registration)
                events_b, movement_b = load_two_env_session(
                    env_b_path, session_ind=2*(day),
                    cell_registration=cell_registration)

                flat_events_a, bins_a, velocity_a = flat_session(
                    events_a, movement_a)
                flat_events_b, bins_b, velocity_b = flat_session(
                    events_b, movement_b)

                # Find active cells
                active_cells_a = np.reshape(
                    np.argwhere(np.sum(flat_events_a, axis=1) > 0), -1)
                active_cells_b = np.reshape(
                    np.argwhere(np.sum(flat_events_b, axis=1) > 0), -1)

                # Find the edge cells of the session
                edge_cells_a = np.reshape(
                    find_env_edge_cells(events_a, movement_a,
                                        min_number_of_events=min_number_of_events), -1,
                    )
                edge_cells_b = np.reshape(
                    find_env_edge_cells(events_b, movement_b,
                                        min_number_of_events=min_number_of_events), -1,
                    )

                # Define the place cells of the session
                place_cells_a = env_a_place_cells[day]
                place_cells_b = env_b_place_cells[day]

                # Find edge place cells
                velocity_indices_a = np.abs(velocity_a) >= velocity_threshold
                velocity_indices_b = np.abs(velocity_b) >= velocity_threshold
                edge_place_cells_a = find_edge_place_cells(
                    flat_events_a[:, velocity_indices_a],
                    bins_a[velocity_indices_a], EDGE_BINS, place_cells_a)
                edge_place_cells_b = find_edge_place_cells(
                    flat_events_b[:, velocity_indices_b],
                    bins_b[velocity_indices_b], EDGE_BINS, place_cells_b)

                # for cell in edge_cells_a:
                #     plot_cell_activity_in_session(
                #         flat_events_a[cell, :], bins_a, "cell no. %d" %cell,
                #         velocity=velocity_a)
                cell_groups = \
                    [active_cells_a, edge_cells_a, place_cells_a,
                     edge_place_cells_a,
                     active_cells_b, edge_cells_b, place_cells_b,
                     edge_place_cells_b]

                cells_overlap = population.iterate_over_sessions(cell_groups,
                                                                 tools.intersection_overlap)
                cells_recurrence = population.iterate_over_sessions(cell_groups,
                                                                    tools.intersection_recurrence)

                plot_matrix(cells_overlap, 'Overlap between groups of cells',
                            "%s day %d" % (mouse, day))
                plt.savefig(os.path.join(OUTPUT_PATH,
                                         'overlap %s day %d.png' % (
                                         mouse, day)))
                plot_matrix(cells_recurrence,
                            'Recurrence between groups of cells',
                            "%s day %d" % (mouse, day))
                plt.savefig(os.path.join(OUTPUT_PATH,
                                         'recurrence %s day %d.png' % (
                                         mouse, day)))
                all_days_overlap.append(cells_overlap)
                all_days_recurrence.append(cells_recurrence)

            except Exception as e:
                print "Error in day %d" % day
                print e

        mean_overlap = np.mean(np.stack(all_days_overlap), axis=0)
        plot_matrix(mean_overlap, 'Overlap between groups of cells', mouse)
        plt.savefig(
            os.path.join(
                OUTPUT_PATH, 'all days mean overlap %s.png' % (mouse)))
        mean_recurrence = np.mean(np.stack(all_days_recurrence), axis=0)
        plot_matrix(mean_recurrence, 'Recurrnce between groups of cells', mouse)
        plt.savefig(
            os.path.join(OUTPUT_PATH,
                         'all days mean recurrence %s.png' % (mouse)))

    raw_input('press enter')

if __name__ == '__main__':
    main()