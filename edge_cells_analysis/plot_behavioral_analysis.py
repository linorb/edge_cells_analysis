import csv
import os
import matplotlib.pyplot as plt
import numpy as np

analysis_path = r"Z:\experiments\projects\positive_reward\analysis"
citric_acid_path = os.path.join(analysis_path, "citric_acid_group")
water_deprived_path = os.path.join(analysis_path, "water_deprived_group")
ensure_path = os.path.join(analysis_path, "ensure_group")

example_day = "day_3"


def load_tracking_file(filename):
    full_tracking = csv.DictReader(open(filename))
    bins = []
    velocity = []
    x = []
    for row in full_tracking:
        bins.append(int(float(row['bin_x'])))
        velocity.append(float(row['velocity_x']))
        x.append(float(row['normalized_centroid_x']))

    tracking = {'bin': np.array(bins), 'velocity': np.array(velocity),
                'x': np.array(x)}

    return tracking


def plot_mice_box_plot(data, field_name, group_name):
    # Data is a list at length of number of mice of the same field (example:
    # a list of velocities across a session
    f, axx = plt.subplots(1, 1)
    axx.boxplot(data)
    axx.set_xlabel("Mice", fontsize=18)
    axx.set_ylabel(field_name, fontsize=18)
    axx.set_title(group_name, fontsize=18)
    axx.set_ylim(-0.1, 7)
    f.show()


def plot_behavioral_trace(location, title):
    Fs = 1/30.0
    time = Fs*np.arange(len(location))
    f, axx = plt.subplots(1, 1)
    axx.plot(time, location)
    axx.set_xlabel("Time [sec]", fontsize=18)
    axx.set_ylabel("Location [bins]", fontsize=18)
    axx.set_ylim(-0.1, 9.1)
    axx.set_title(title, fontsize=18)
    f.show()


def load_group_behavior(group_analysis_path):
    mice_names = os.listdir(group_analysis_path)
    mice_behavior = []
    for mouse in mice_names:
        filename = os.path.join(group_analysis_path, mouse, example_day,
                                'tracking', 'behavior.csv')
        mice_behavior.append(load_tracking_file(filename))

    return mice_behavior


def create_mice_speed_list(mice_behavior):
    # Calculate the normalized velocity and multiply by 60 (cm)
    mice_field = []
    for mouse in mice_behavior:
        location = mouse['bin'][:-1]
        x = mouse['x']
        speed = np.abs((x[1:] - x[:-1])*60)
        # Taking only the velocity *not* on edges
        mice_field.append(speed[(location > 0) & (location < 9)])
    return mice_field

def main():
    ca_mice = load_group_behavior(citric_acid_path)
    wd_mice = load_group_behavior(water_deprived_path)
    en_mice = load_group_behavior(ensure_path)

    ca_speed = create_mice_speed_list(ca_mice)
    wd_speed = create_mice_speed_list(wd_mice)
    en_speed = create_mice_speed_list(en_mice)

    plot_mice_box_plot(ca_speed, 'Speed', 'Citric acid group')
    plot_mice_box_plot(wd_speed, 'Speed', 'Water deprived group')
    plot_mice_box_plot(en_speed, 'Speed', 'Ensure group')

    plot_behavioral_trace(ca_mice[1]['bin'], 'Citric acid example mouse')
    plot_behavioral_trace(wd_mice[4]['bin'], 'Water deprived example mouse')
    plot_behavioral_trace(en_mice[3]['bin'], 'Ensure example mouse')
    
    raw_input('press enter')

if __name__ == '__main__':
    main()


