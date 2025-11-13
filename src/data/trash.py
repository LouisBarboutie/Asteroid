from matplotlib.figure import Figure
import time



# Load the data from the ephemeris
kernel = SPK.open('data/de421.bsp')
# Make vector output attractive
np.set_printoptions(precision=3)


# Position in km, velocity in km/day
mars_position, mars_velocity = kernel[0,4].compute_and_differentiate(date)
mars_velocity_seconds = mars_velocity/86400.0

#print(mars_position, mars_velocity, mars_velocity_seconds)

kernel.close()
try:
    element_data = data['result']['table'][0]
    eccentricity = element_data['EC']
    semimajor_axis = element_data['AC']
    inclination = element_data['IN']
    print(f"Success")
    print(f"Eccentricity: {eccentricity}")
    print(f"Semimajor-axis: {semimajor_axis}")
    print(f"Inclination: {inclination}")
except KeyError as e:
    print(f"error extracting data")
    if 'error' in data:
        print(f"API error message: {data['error']}")
try:
    # Get the controller for 'firefox'
    browser = webbrowser.get('firefox')

except webbrowser.Error:
    print(f"Could not find Firefox. Plot saved to {file_path}")
    # Fallback to default browser if Firefox isn't found
    # webbrowser.open('file://' + os.path.realpath(file_path))