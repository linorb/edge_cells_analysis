"""This analysis looks at the recurrence of place cells\edge cells in
one environment vs. another"""
from cell_activity_in_different_environments import *
from zivlab.analysis.place_cells import find_place_cells
from zivlab.analysis import population
from bambi_playback_analysis.code import tools

VELOCITY_THRESHOLD = 5


def find_edge_place_cells(events, bins, edge_bins, place_cells_indices):
    """Finds the place cells that represents the edge out of the place cells"""
    event_rate_distribution = population._calculate_event_rate_distribution(
        bins, events[place_cells_indices, :])
    bin_tuning = np.argmax(event_rate_distribution, axis=1)
    edge_cell_indices = []
    for edge in edge_bins:
        edge_cell_indices.extend(np.argwhere(bin_tuning == edge))
    edge_cell_indices = np.reshape(np.array(edge_cell_indices), -1)
    return place_cells_indices[edge_cell_indices]


def plot_overlap_matrix(overlap, title):
    f, axx = plt.subplots(1, 1)
    axx.imshow(overlap, aspect='auto', interpolation='none')
    axx.set_xticks(np.arange(0, 8, 1))
    axx.set_yticks(np.arange(0, 8, 1))
    axx.set_xticklabels(['AC_a', 'EC_a', 'PC_a', 'EPC_a',
                         'AC_b', 'EC_b', 'PC_b', 'EPC_b'])
    axx.set_yticklabels(['AC_a', 'EC_a', 'PC_a', 'EPC_a',
                         'AC_b', 'EC_b', 'PC_b', 'EPC_b'])
    axx.set_title('Overlap between group of cells', fontsize=18)
    plt.text(0, 8.5, "AC = Active Cells, EC = Edge Cells, PC = Place cells, "
                     "EPC = Edge Place Cells", size=10)
    plt.draw()
    f.colorbar
    f.suptitle(title, fontsize=18)
    f.show()
    
    return


def main():
    mice_data = dict.fromkeys(MICE)
    for mouse in MICE:
        print mouse
        all_days_overlap = []
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
                active_cells_a = np.reshape(np.argwhere(np.sum(flat_events_a, axis=1) > 0), -1)
                active_cells_b = np.reshape(np.argwhere(np.sum(flat_events_b, axis=1) > 0), -1)

                # Find the edge cells of the session
                edge_cells_a = np.reshape(find_env_edge_cells(events_a, movement_a), -1)
                edge_cells_b = np.reshape(find_env_edge_cells(events_b, movement_b), -1)

                # Find place cells
                # TODO: find out the reason for the large number of place cells - maybe take out the edges?
                velocity_indices_a = np.abs(velocity_a) >= VELOCITY_THRESHOLD
                place_cells_a, _, _ = find_place_cells(
                    bins_a[velocity_indices_a], flat_events_a[:, velocity_indices_a],
                    min_number_of_events=5)

                velocity_indices_b = np.abs(velocity_b) >= VELOCITY_THRESHOLD
                place_cells_b, _, _ = find_place_cells(
                    bins_b[velocity_indices_b], flat_events_b[:, velocity_indices_b],
                    min_number_of_events=5)

                # Find edge place cells
                edge_place_cells_a = find_edge_place_cells(
                    flat_events_a[:, velocity_indices_a],
                    bins_a[velocity_indices_a], EDGE_BINS, place_cells_a)
                edge_place_cells_b = find_edge_place_cells(
                    flat_events_b[:, velocity_indices_b],
                    bins_b[velocity_indices_b], EDGE_BINS, place_cells_b)

                cell_groups = \
                    [active_cells_a, edge_cells_a, place_cells_a, edge_place_cells_a,
                     active_cells_b, edge_cells_b, place_cells_b, edge_place_cells_b]
                
                cells_overlap = population.iterate_over_sessions(cell_groups,
                    tools.intersection_overlap)
                
                import pdb; pdb.set_trace()
                plot_overlap_matrix(cells_overlap, "%s day %d" % (mouse, day))
                all_days_overlap.append(cells_overlap)

            except Exception as e:
                print "Error in day %d" % day
                print e
                import pdb; pdb.set_trace()
        
        import pdb; pdb.set_trace()
        mean_overlap = np.mean(np.stack(all_days_overlap), axis=2)
        plot_overlap_matrix(mean_overlap, mouse)


if __name__ == '__main__':
    main()