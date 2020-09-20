import math
import numpy as np
import eel

eel.init('web')


@eel.expose
def py_func(str):
    str = str.split(' ')
    str = [float(i) for i in str]
    id, e, t_oa, i_k, omega_point, sqrt_A, omega_0, w, M_0 = str
    # id = 8
    # e = 0.4450321198E-002  # Eccentricity
    # t_oa = 405504.0  # Time of Applicability(sec)
    # i_k = 0.9705328666  # Orbital Inclination(rad)
    # omega_point = -0.7886042771E-008  # Rate of Right Ascen(r/s)
    # sqrt_A = 5153.510742  # SQRT(A)  (m**(1/2)) - корень из значения большой полуоси
    # omega_0 = 0.2309270258E+000  # Right Ascen at Week(rad)
    # w = -0.232577104  # Argument of Perigee(rad)
    # M_0 = -0.4053031157E+000  # Mean Anom(rad)

    # Constants
    mju = 3.986005 * 10 ** 14  # WGS 84 value of the earth's gravitational constant for GPS user
    omega_e_point = 7.2921151467 * (10 ** (-5))  # WGS 84 value of the earth's rotation rate

    # Elements of Coordinate Systems (IS-GPS-200F.pdf sheet 105)
    t0 = 3 * 24 * 3600 + id * 3600  # Current time (8:00 am, 18.09.2019) math.since beginning of the week (sec) = 3 days * 24 hours * 3600 sec + 8 * 3600 sec
    dt = 18  # Deviation current time from GPS time
    t = t0 + dt  # Current time GPS (sec)
    t_k = t - t_oa  # Time from ephemeris reference epoch (sec)
    toe = 604800  # The epoch time

    # Condition:
    # tk shall be the actual total time difference between the time t and the...
    # epoch time toe, and must account for beginning or end of week crossovers.
    # That is, if tk is greater than 302,400 seconds, subtract 604,800 seconds...
    # from tk. If tk is less than -302,400 seconds, add 604,800 seconds to tk.

    if t_k < -302400:
        t_k = t_k + toe
    elif t_k > 302400:
        t_k = t_k - toe

    A = sqrt_A ** 2  # Semi-major axis (m)
    n_0 = math.sqrt(mju / A ** 3)  # Computed mean motion (rad/sec)
    n = n_0  # Corrected mean motion
    M_k = M_0 + n * t_k  # Mean anomaly (rad

    # Finding Ek

    f = lambda Ek1, Mk, e: Mk - (Ek1 - e * math.sin(Ek1))  # Kepler's equation Mk = Ek - e*math.sin(Ek)

    Ek_0 = M_k - e  # minimum of Ek
    eps = 0.0000000001  # required accuracy
    NN = 1001  # number of iterations

    def g(Ek1, Mk, e):
        return Mk + e * math.sin(Ek1)

    E_k = g(Ek_0, M_k, e)

    for i in range(1, NN):
        if abs(E_k - Ek_0) <= eps:
            break
        Ek_0 = E_k
        E_k = g(Ek_0, M_k, e)  # Eccentric Anomaly

    nju_k = math.atan2((math.sqrt(1 - e ** 2) * math.sin(E_k) / (1 - e * math.cos(E_k))),
                       ((math.cos(E_k) - e) / (1 - e * math.cos(E_k))))  # True anomaly

    F_k = nju_k + w  # Argument of Latitude (rad)
    u_k = F_k  # Corrected Argument of Latitude (rad)
    r_k = A * (1 - e * math.cos(E_k))  # Corrected Radius (rad)

    # Positions in orbital plane:
    x_k_hatch = r_k * math.cos(u_k)  # (m)
    y_k_hatch = r_k * math.sin(u_k)  # (m)

    omega_k = omega_0 + (
            omega_point - omega_e_point) * t_k - omega_e_point * t_oa  # Corrected longitude of ascending node (rad)

    # Earth Fixed coordinates:
    x_k = x_k_hatch * math.cos(omega_k) - y_k_hatch * math.cos(i_k) * math.sin(omega_k)  # (m)
    y_k = x_k_hatch * math.sin(omega_k) + y_k_hatch * math.cos(i_k) * math.cos(omega_k)  # (m)
    z_k = y_k_hatch * math.sin(i_k)  # (m)

    # Bauman University roof coordinates:
    x = 2846228.896  # (m)
    y = 2198658.103  # (m)
    z = 5249983.343  # (m)

    # Counting athimuth and elevation
    # The coordinates of the navigation spacecraft range vector with respect to the roof:
    delta_r_x = x_k - x  # (m)
    delta_r_y = y_k - y  # (m)
    delta_r_z = z_k - z  # (m)
    delta_r = math.sqrt(delta_r_x ** 2 + delta_r_y ** 2 + delta_r_z ** 2)  # range vector (m)

    a = 6378137  # semi-major axis of an ellipsoid
    b = 6356752.3142  # minor axis of an ellipsoid

    LL = math.atan2(y, x)  # longotude
    LL_grad = LL * 180 / math.pi

    # Medvedev algorithm:
    ff = 298.257223  # compression
    k0 = (ff - 1) / ff
    k1 = a * (2 * ff - 1) / ff / (ff - 1)
    k2 = k0 * k1
    R = math.sqrt(x ** 2 + y ** 2)
    U = math.atan((k1 / math.sqrt(z ** 2 + (k0 * R) ** 2) + 1) * k0 * z / R)
    BB = math.atan2(z + k1 * (math.sin(U)) ** 3, R - k2 * (math.cos(U)) ** 3)  # latitude
    BB_grad = BB * 180 / math.pi
    H = R * math.cos(BB) + z * math.sin(BB) - a * math.sqrt(1 - k2 * (math.sin(BB)) ** 2 / a)  # latitude

    M1 = np.array([[delta_r_x], [delta_r_y], [delta_r_z]])
    # Move on topografic rectangular coorsinates
    M2 = np.array([[-math.sin(BB) * math.cos(LL), -math.sin(BB) * math.sin(LL), math.cos(BB)],
                   [-math.sin(LL), math.cos(LL), 0.],
                   [math.cos(BB) * math.cos(LL), math.cos(BB) * math.sin(LL), math.sin(BB)]])

    M3 = M2.dot(M1)

    # zenith distance:

    Z_r = math.atan2(np.sqrt(M3[0] ** 2 + M3[1] ** 2), M3[2])
    Z_r_grad = Z_r * 180 / math.pi

    beta = math.pi / 2 - Z_r  # evidence (rad)
    alpha = math.atan2(M3[1], M3[0])  # azimuth (rad)

    beta_grad = beta * 180 / math.pi  # evidence (grad)
    # print(beta)
    alpha_grad = alpha * 180 / math.pi  # azimuth (grad)

    # Notice that azimuth is in range (0,360) grad

    if alpha_grad < 0:
        alpha_grad = alpha_grad + 360
    elif alpha_grad > 360:
        alpha_grad = alpha_grad - 360

    # print(f'Satelite coordinates: {x_k,y_k,z_k} \n')
    # print(f'Distance: {math.sqrt(x_k**2+y_k**2+z_k**2)} m\n')
    # print(f'Azimuth: { alpha_grad} grad\n')
    # print(f'Elevation: {beta_grad} grad\n')
    # print('---------------\n')
    return f'{x_k} {y_k} {z_k} {math.sqrt(x_k ** 2 + y_k ** 2 + z_k ** 2)} {alpha_grad} {beta_grad}'


eel.start('main2.html', size=(1000, 700))
