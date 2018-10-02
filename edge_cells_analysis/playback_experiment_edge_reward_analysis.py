"""This code look at the rewarded and not rewarded epochs at the playback
experiment in different analysis methods"""
import numpy as np
import os
import matplotlib.pyplot as plt
from cPickle import load
from bambi_playback_analysis.code.dynamics import load_cell_registration
from bambi_playback_analysis.code import tools
from bambi.tools.activity_loading import unite_sessions
from zivlab.analysis import population


BAMBI_MICE = ['c51m4', 'c58m5', 'c59m1', 'c59m9', 'c60m5']
SESSIONS = np.arange(6)
ANALYSIS_PATH = R'Z:\experiments\projects\bambi\linear_track_2\analysis'
EXPERIMENT_PHASE = '2_groups'
EDGE_BINS = [0, 9]


def load_session(mouse, session):
    cell_reg_filename = os.path.join(ANALYSIS_PATH, 'registration',
                                     'sessions_6-12', mouse,
                                     'cellRegistered.mat')
    cell_registration = load_cell_registration(cell_reg_filename)

    phase_path = os.path.join(ANALYSIS_PATH, EXPERIMENT_PHASE)
    session_path = os.path.join(phase_path, mouse, 'day_%d' % (session+1),
                                'cooked')

    events = tools.load_events(
        os.path.join(session_path, 'finalEventsMat.mat'))
    behavior = tools.load_tracking_file(
        (os.path.join(session_path, 'behavior.csv')))
    global_events_numbering = unite_sessions(
        [events], [session], cell_registration)

    water_trace = load(open(os.path.join(session_path,
                                                 'water_system_trace.pkl')))
    water_reward_frames = np.array([x[2] for x in water_trace[1] if x[2]>0])

    return global_events_numbering, behavior, water_reward_frames


def segment_run_edge(bins, edge_bins):
    """Segment the bins to edge and run bins. return the number of frames for
    each segment"""
    edge_bins_mask = np.zeros_like(bins, dtype=bool)
    for b in edge_bins:
        edge_bins_mask[bins == b] = True

    edge_masked = np.ma.array(bins, mask=edge_bins_mask)
    run_locations = np.ma.flatnotmasked_contiguous(edge_masked)

    run_masked = np.ma.array(bins, mask=~edge_bins_mask)
    edge_locations = np.ma.flatnotmasked_contiguous(run_masked)

    return run_locations, edge_locations


def tag_reward(segment_indices, reward_indices):
    """Marks the segments according to the reward dispensing indices. A segment
    which contains a dispensing frame will be marked as rewarded - 1,
    returns segments reward"""
    reward_segments = []
    for reward in reward_indices:
        for segment in segment_indices:
                if (reward >= segment.start) & (reward < segment.stop):
                    reward_segments.append(segment)
                    break

    not_reward_segments = []
    for segment in segment_indices:
        if not(np.isin(segment, reward_segments)):
            not_reward_segments.append(segment)

    return reward_segments, not_reward_segments


def main():
    for mouse in BAMBI_MICE:
        for session_ind in SESSIONS:
            # Load session
            events , behavior, water_reward_frames = load_session(mouse, session_ind)
            bins = behavior['bin']

            # Create segments of (not) rewarded edge phases
            _, edge_locations = segment_run_edge(bins, EDGE_BINS)
            reward_segments, not_reward_segments = tag_reward(
                edge_locations, water_reward_frames)
            reward_events = [events[:, segment] for segment in reward_segments]
            non_reward_events = [events[:, segment] for segment in not_reward_segments]
            number_of_rewarded_epochs = len(reward_segments)

            events_segments = reward_events
            events_segments.extend(non_reward_events)

            # Calculate the correlation between types of edge phases
            ensemble_correlation_matrix = population.iterate_over_sessions(
                events_segments, population.ensemble_correlation)

            f, axx = plt.subplots(1, 1)
            im = axx.imshow(ensemble_correlation_matrix, aspect='auto',
                            interpolation='none')
            axx.axhline(y=number_of_rewarded_epochs, linewidth=3, color='k')
            axx.axvline(x=number_of_rewarded_epochs, linewidth=3, color='k')
            axx.set_title('%s day %d' % (mouse, session_ind + 1), fontsize=18)
            f.suptitle('Ensemble correlation between edge epochs', fontsize=18)
            f.show()


if __name__ == '__main__':
    main()

