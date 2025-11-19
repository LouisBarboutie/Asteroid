def angular_momentum_normalized(i, raan):
    h = [np.sin(i)*np.sin(raan), -np.sin(i)*np.cos(raan), np.cos(i)]
    return h
def nu_to_eccentric_anomaly(e, nu):
    E = 2*np.arctan(np.tan(nu/2)*np.sqrt((1-e)/(1+e)))
    return E
def time_to_raan(w, nu, mu, a, e):
    # nu corresponds to current true anomaly
    E0 = nu_to_eccentric_anomaly(e, nu) # Initial E
    E = nu_to_eccentric_anomaly(e, -w) # Final E
    dt = np.sqrt(a**3/mu)*((E-e*np.sin(E))-(E0-e*np.sin(E0)))

    if dt <= 0:
        dt = dt + 2*np.pi*np.sqrt(a**3/mu)
        return dt
    else:
        return dt


h_earth = angular_momentum_normalized(elements_earth[2], elements_earth[3])
h_asteroid = angular_momentum_normalized(elements_asteroid[2], elements_asteroid[3])
relative_raan_vector = np.cross(h_earth, h_asteroid)
relative_raan_vector_normalized = relative_raan_vector / np.sqrt(
relative_raan_vector[0] ** 2 + relative_raan_vector[1] ** 2 + relative_raan_vector[2] ** 2)
relative_raan_angle = np.arctan2(relative_raan_vector[0], relative_raan_vector[1])

dt = time_to_raan(elements_earth[4], elements_earth[5], mu_sun, elements_earth[0], elements_earth[1])
