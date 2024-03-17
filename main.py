import csv
import matplotlib.pyplot as plt
import mpl_axes_aligner

class LinearInterpolator:
    def __init__(self, csv_file):
        self.xs, self.ys = self.read_csv(csv_file)

    def read_csv(self, csv_file):
        xs = []
        ys = []
        
        with open(csv_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                xs.append(float(row[0]))
                ys.append(float(row[1]))
                
        return xs, ys

    def interpolate(self, x):
        if x > max(self.xs):
            return 0
        
        for idx, lx in enumerate(self.xs):
            if x == lx:
                return self.ys[idx]

            elif x < lx:
                break

        idx_up = idx
        idx_dn = idx - 1

        x1 = self.xs[idx_dn]
        x2 = self.xs[idx_up]

        y1 = self.ys[idx_dn]
        y2 = self.ys[idx_up]
        # Linear interpolation formula
        y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)

        return y

    def max_x(self):
        return max(self.xs)

def sign(x): return 1 if x >= 0 else -1

def plot_traj(xs, hs, vels, accels, times):
    fig, ax = plt.subplots()

    plt.title("Trajectory")
    ax.plot(times, hs, label="Alt. AGL")
    ax.plot(times, vels, label="Vel.")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Alt. AGL (m), Vel. (m/s)")
    plt.grid()
    ax2 = ax.twinx()
    ax2.plot(times, accels, label="Accel.", color="green")
    ax2.set_ylabel("Accel. (m/s)")
    mpl_axes_aligner.align.yaxes(ax, 0, ax2, 0, 0.2)

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)
    
    plt.show()

def plot_forces(times, thrusts, drags, gravs):
    plt.figure(2)
    plt.title("Engine Thrust Profile")
    plt.plot(times, thrusts)
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.grid()

    plt.figure(3)
    plt.title("Drag Profile")
    plt.plot(times, drags)
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.grid()

    plt.figure(4)
    plt.title("Gravitational Acceleration Profile")
    plt.plot(times, gravs)
    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (m s-2)")
    plt.grid()

    plt.show()

def main(x0, drag_coeff, cross_sec, mass_init, mass_final):
    density_interpolator = LinearInterpolator('./data/atm_density.csv')
    thrust_interpolator = LinearInterpolator('./data/A8-3.csv')
    burn_time = thrust_interpolator.max_x()
    mass_flow_avg = (mass_init - mass_final) / burn_time

    def alt2grav(altitude):
        gravity = 9.80665 * (6369000/(6369000+altitude))**2
        return gravity

    def calc_drag(velocity, altitude):
        density = density_interpolator.interpolate(altitude)
        drag = (0.5 * density * velocity**2 * drag_coeff * cross_sec) * -sign(velocity)
        return drag

    xs = []
    vels = []
    accels = []
    hs = []
    times = []
    thrusts = []
    drags = []
    gravs = []

    x = x0 # altitude ASL (m)
    v = 0 # velocity (m s-1)
    a = 0 # acceleration (m s-2)

    m = mass_init # mass (kg)

    h = 0 # altitude AGL (m)
    t = 0 # time (s)

    dt = 0.001

    cycle = 0
    running = True
    stage = "LAUNCH" # launch, coast, landing
    while running:

        if stage == "LAUNCH":
            thrust = thrust_interpolator.interpolate(t)
            drag = calc_drag(v, x) # will come out negative when going up
            grav = alt2grav(x)
            a = (thrust + drag) / m - grav

            v += a * dt
            x += v * dt
            h = x - x0

            m = m - mass_flow_avg * dt

            # not enough thrust to launch yet
            if x <= x0 and v <= 0:
                x = x0
                v = 0
                a = 0
                h = 0

            if t > burn_time:
                stage = "COAST"

        elif stage == "COAST":
            drag = calc_drag(v, x)
            grav = alt2grav(x)
            a = drag / m - grav

            v += a * dt
            x += v * dt
            h = x - x0

        # record data
        xs.append(x)
        vels.append(v)
        accels.append(a)
        hs.append(h)
        times.append(t)
        thrusts.append(thrust)
        drags.append(drag)
        gravs.append(grav)

        # end conditions
        # hit ground
        if stage != "LAUNCH" and x <= x0:
            running = False
        
        cycle += 1
        t = cycle * dt

    plot_traj(xs, hs, vels, accels, times)
    plot_forces(times, thrusts, drags, gravs)

    return xs, hs, vels, accels, times

def init():
    ## ENTER INPUTS HERE
    x0 = # m
    drag_coeff =
    cross_sec =  # m2
    mass_init =  # kg
    mass_final =  # kg
    main(x0, drag_coeff, cross_sec, mass_init, mass_final)

if __name__ == "__main__":
    init()
