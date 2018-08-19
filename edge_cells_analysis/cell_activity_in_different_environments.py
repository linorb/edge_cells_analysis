""" This analysis looks at the cell level on the activity of edge cells
on different tracks (L-shape and linear)
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from bambi_playback_analysis.code.long_linear_edge_cells import find_edge_cells
from bambi.tools import matlab, activity_loading

EDGE_BINS = [0, 1, 10, 11]
MIN_NUMBER_OF_EVENTS = 10
EDGE_THRESHOLD = 0.5
FRAME_RATE = 20
DATA_PATH = R'D:\dev\replays\work_data\two_environments'
MICE = ['c6m4', 'c7m4', 'c11m1', 'c13m1']
DAYS = np.arange(1, 8)
CELL_REGISTRATION_PATH = \
    R"D:\dev\replays\work_data\two_environments"


def load_two_env_session(session_path, session_ind=np.nan, cell_registration=[]):
    """Load session. output:
    linear_trials_events: A list of events matrices, where each element contains
        the events matrix of a single trial (ordered by their index)
    movement_data: A list of the movment information (filname,x,y,bin,velocity)
        from all trials"""
    # Load events (the fixed matrix marks the events from the rise time)
    events_filename = 'FixedEventsMat.mat'
    log_filename = 'frameLog.csv'
    all_events = matlab.load_events_file(
        os.path.join(session_path, events_filename)) > 0
    frame_log = matlab.load_frame_log_file(
        os.path.join(session_path, log_filename))
    # Register cells
    if not(cell_registration==[]):
        all_events = activity_loading.unite_sessions(
            [all_events], [session_ind], cell_registration)
    # Take only the linear track events
    linear_trials_events = activity_loading.order_events_into_trials(
        all_events, frame_log[1:-1])

    # Load behavior
    if 'envB' in session_path:
        behavior_filename = 'my_mvmt_fixed.mat'
    else:  # envA or linear
        behavior_filename = 'my_mvmt_smooth.mat'
    movement_data = matlab.load_mvmt_file(
        os.path.join(session_path, behavior_filename))[1:]
    # widening the bins - to create a total of 12 bins instead of 24
    for trial_behavior in movement_data:
        wider_bins = activity_loading.wide_binning(
            trial_behavior['bin'], 24, 2)
        trial_behavior['bin'] = wider_bins

    return linear_trials_events, movement_data


def flat_session(events, movement_data):
    """Flattens the trials division of the bins, velocity and events of the
    session. Returns the events as one matrix and the bins & velocity as one
    array each"""
    flat_events = np.hstack(events)
    bins = []
    velocity = []
    for trial in movement_data:
        bins.append(trial['bin'][1:])
        velocity.append(trial['velocity'][1:])
    flat_bins = np.hstack(bins)
    flat_velocity = np.hstack(velocity)

    return flat_events, flat_bins, flat_velocity


def plot_comparison_of_two_sessions(cell_events_a, bins_a,
                                    cell_events_b, bins_b, title):
    """Plot a comparison of events between two sessions of two session from
    different environments"""
    f, axx = plt.subplots(1, 2)
    axx[0].plot(bins_a, 'b')
    cell_events_a[cell_events_a == 0] = np.nan
    axx[0].plot(bins_a*cell_events_a, '*r', markersize=12)
    axx[0].set_xlabel('Time', fontsize=18)
    axx[0].set_ylabel('Bins', fontsize=18)
    axx[0].set_title('Environment A events', fontsize=18)

    axx[1].plot(bins_b, 'b')
    cell_events_b[cell_events_b == 0] = np.nan
    axx[1].plot(bins_b * cell_events_b, '*r', markersize=12)
    axx[1].set_xlabel('Time', fontsize=18)
    axx[1].set_ylabel('Bins', fontsize=18)
    axx[1].set_title('Environment B events', fontsize=18)

    f.suptitle(title, fontsize=20)
    f.show()

    return


def plot_hist_comparison(cell_events_a, bins_a,
                         cell_events_b, bins_b, title):
    number_of_bins = np.max(bins_a)
    f, axx = plt.subplots(1, 2)
    cell_events_a[cell_events_a == 0] = np.nan
    events_count_a = bins_a * cell_events_a
    n_a, _, _ = axx[0].hist(events_count_a[~np.isnan(events_count_a)],
                bins=np.arange(-0.1, number_of_bins+0.1))
    axx[0].set_ylabel('Event count', fontsize=18)
    axx[0].set_xlabel('Bins', fontsize=18)
    axx[0].set_xlim(-0.1, number_of_bins+1)
    axx[0].set_title('Environment A events', fontsize=18)

    cell_events_b[cell_events_b == 0] = np.nan
    events_count_b = bins_b * cell_events_b
    n_b, _, _ = axx[1].hist(events_count_b[~np.isnan(events_count_b)],
                bins=np.arange(-0.1, number_of_bins+0.1))
    axx[1].set_xlabel('Bins', fontsize=18)
    axx[1].set_xlim(-0.1, number_of_bins+1)
    axx[1].set_title('Environment B events', fontsize=18)

    f.suptitle(title, fontsize=20)
    # f.show()

    return n_a, n_b


def plot_tuning_curve_comparison(cell_events_a, bins_a, cell_events_b,
                                 bins_b, title):
    bins_range = np.arange(-0.5, np.max(bins_a)+1.5)
    f, axx = plt.subplots(1, 2)
    cell_events_a[cell_events_a == 0] = np.nan
    where_events = bins_a * cell_events_a
    events_count_a, _ = np.histogram(where_events[~np.isnan(where_events)],
                                     bins=bins_range)
    occupancy_a, _ = np.array(np.histogram(bins_a, bins=bins_range))
    occupancy_a = np.array(occupancy_a, dtype=float)
    tuning_curve_a = events_count_a/occupancy_a
    axx[0].bar(bins_range[:-1], tuning_curve_a)
    axx[0].set_ylabel('P(event|bin)')
    axx[0].set_xlabel('Bins', fontsize=18)
    axx[0].set_xlim(-0.5, np.max(bins_a)+1)
    axx[0].set_title('Environment A events', fontsize=18)

    cell_events_b[cell_events_b == 0] = np.nan
    where_events = bins_b * cell_events_b
    events_count_b, _ = np.histogram(where_events[~np.isnan(where_events)],
                                     bins=bins_range)
    occupancy_b, _ = np.histogram(bins_b, bins=bins_range)
    occupancy_b = np.array(occupancy_b, dtype=float)
    tuning_curve_b = events_count_b / occupancy_b
    axx[1].bar(bins_range[:-1], tuning_curve_b)
    axx[1].set_ylabel('P(event|bin)')
    axx[1].set_xlabel('Bins', fontsize=18)
    axx[1].set_xlim(-0.5, np.max(bins_b) + 1)
    axx[1].set_title('Environment B events', fontsize=18)

    f.suptitle(title, fontsize=20)
    f.show()

    return tuning_curve_a, tuning_curve_b


def find_env_edge_cells(events, movement):
    trial_bins = []
    for trial in movement:
        trial_bins.append(trial['bin'][1:])

    return find_edge_cells(
        events, trial_bins, MIN_NUMBER_OF_EVENTS, EDGE_THRESHOLD,
        edge_bins=EDGE_BINS, frame_rate=FRAME_RATE)


def plot_mouse_summary(number_of_edge_cells_a, number_of_edge_cells_b, 
                       unique_edge_cells, edge_cells_corr, mouse_name):
    f, axx = plt.subplots(1, 2)
    axx[0].scatter(number_of_edge_cells_a, number_of_edge_cells_b)
    axx[0].set_xlabel('number of edge cells environment a', fontsize=18)
    axx[0].set_ylabel('number of edge cells environment b', fontsize=18)
    axx[0].set_title('Number of edge cells', fontsize=18)
    
    axx[1].scatter(unique_edge_cells, edge_cells_corr)
    axx[1].set_xlabel('Number of unique edge cells', fontsize=18)
    axx[1].set_ylabel('Fraction of high correlation (>0.65) cells', fontsize=18)
    axx[1].set_title('High correlation edge cells', fontsize=18)
    
    f.suptitle(mouse_name, fontsize=20)
    f.show()
    
    return


def main():
    for mouse in MICE:
        print mouse
        number_of_edge_cells_a = []
        number_of_edge_cells_b = []
        unique_edge_cells = []
        edge_cells_corr = []
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
                    env_a_path, session_ind=2*(day-1), cell_registration=cell_registration)
                events_b, movement_b = load_two_env_session(
                    env_b_path, session_ind=2*(day), cell_registration=cell_registration)

                # Find the edge cells of the session
                edge_cells_a = find_env_edge_cells(events_a, movement_a)
                edge_cells_b = find_env_edge_cells(events_b, movement_b)
                all_edge_cells = np.unique(np.concatenate([edge_cells_a, edge_cells_b]))

                # Plot edge cells activity
                flat_events_a, flat_bins_a, _ = flat_session(events_a, movement_a)
                flat_events_b, flat_bins_b, _ = flat_session(events_b, movement_b)
                cells_event_count = []
                cells_corr = []
                for cell in all_edge_cells:
                    title = "Cell no. %d" % cell
                    cell_activity_a = flat_events_a[cell, :]
                    cell_activity_b = flat_events_b[cell, :]
                    # plot_comparison_of_two_sessions(cell_activity_a, flat_bins_a,
                    #                                 cell_activity_b, flat_bins_b, title)
                    n_a, n_b = \
                        plot_hist_comparison(cell_activity_a, flat_bins_a,
                                                     cell_activity_b, flat_bins_b, title)
                    cells_event_count.append([n_a, n_b])
                    cells_corr.append(np.corrcoef(n_a, n_b)[0, 1])

                    # print "Cell no. %d correlation coefficient: %f" % \
                    #       (cell, cells_corr[-1])
                
                number_of_edge_cells_a.append(len(edge_cells_a))
                number_of_edge_cells_b.append(len(edge_cells_b))
                unique_edge_cells.append(len(all_edge_cells))
                edge_cells_corr.append(np.sum(np.array(
                    cells_corr) > 0.65)/np.float(len(all_edge_cells)))
                print "number of edge cells in environment A: %d" % len(edge_cells_a)
                print "number of edge cells in environment B: %d" % len(edge_cells_b)
                print "number of unique edge cells: %d" % len(all_edge_cells)
                print "Fraction of edge cells above 0.65 correlation: %f" % edge_cells_corr[-1]
                plt.close('all')
            
            except Exception as e:
                print "Error in day %d" % day
                import pdb; pdb.set_trace()
               
        plot_mouse_summary(number_of_edge_cells_a, number_of_edge_cells_b, 
            unique_edge_cells, edge_cells_corr, mouse)
        
        raw_input('press enter')
        
if __name__ == '__main__':
    main()