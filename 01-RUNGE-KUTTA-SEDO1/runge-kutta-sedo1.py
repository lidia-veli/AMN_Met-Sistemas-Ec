def runge_kutta_sedo1(funcX, FuncY, a, b, y0, y1, N):
    '''Metodo Runge Kutta(ord. 4) para sistemas de ecuaciones diferenciales de orden 1'''

    h = (b - a) / N

    t = [a]
    u = [y0]
    v = [y1]

    for k in range(N):
        k11 = funcX(t[k], u[k], v[k])
        k12 = FuncY(t[k], u[k], v[k])

        k21 = funcX(t[k] + h/2, u[k] + h*k11/2, v[k] + h*k12/2)
        k22 = FuncY(t[k] + h/2, u[k] + h*k11/2, v[k] + h*k12/2)

        k31 = funcX(t[k] + h/2, u[k] + h*k21/2, v[k] + h*k22/2)
        k32 = FuncY(t[k] + h/2, u[k] + h*k21/2, v[k] + h*k22/2)

        k41 = funcX(t[k] + h, u[k] + h*k31, v[k] + h*k32)
        k42 = FuncY(t[k] + h, u[k] + h*k31, v[k] + h*k32)

        # añadimos los valores a la lista
        t.append(t[k] + h)
        u.append(u[k] + h*(k11 + 2*k21 + 2*k31 + k41)/6)
        v.append(v[k] + h*(k12 + 2*k22 + 2*k32 + k42)/6)

    return t, u, v
