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
url_mars = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='499'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_mercury = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='199'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_venus = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='299'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)

elements_mercury = get_elements(url_mercury)
elements_venus = get_elements(url_venus)
elements_mars = get_elements(url_mars)
orb_mercury = define_orbit(Sun, elements_mercury)
orb_venus = define_orbit(Sun, elements_venus)

orb_mars = define_orbit(Sun, elements_mars)


#plotter.plot(orb_mercury, label = 'Mercury', color = '#B1B1B1')
#plotter.plot(orb_venus, label = 'Venus', color = '#EEDC82')
#plotter.plot(orb_mars, label = 'Mars', color = '#C1440E')

plotter.plot_ephem(
    orb_earth,
    label = 'Earth',
    color = '#2E8B57'
)

h_earth = angular_momentum_normalized(elements_earth[2],elements_earth[3])
h_asteroid = angular_momentum_normalized(elements_asteroid[2],elements_asteroid[3])
relative_raan_vector = np.cross(h_earth, h_asteroid)

relative_raan_vector_normalized = relative_raan_vector/np.sqrt(relative_raan_vector[0]**2+relative_raan_vector[1]**2+relative_raan_vector[2]**2)
relative_raan_angle = np.arctan2(relative_raan_vector[0], relative_raan_vector[1])*180/np.pi
print(relative_raan_vector_normalized)
print(relative_raan_angle)
print("Line of nodes angle between the two orbits:", relative_raan_angle)
# Determine transfer date for impulsive maneuver
dt = time_to_raan(elements_asteroid_radians[0], elements_asteroid_radians[1], elements_asteroid_radians[4], elements_asteroid_radians[5], mu_sun)



k = np.linspace(-149e6, 149e6, 100)
r_coordinates = np.outer(k, relative_raan_vector_normalized)
r_coordinates_quantity = r_coordinates * u.km
x,y,z = r_coordinates_quantity[0], r_coordinates_quantity[1], r_coordinates_quantity[2]
relative_raan_coordinates = CartesianRepresentation(x,y,z)