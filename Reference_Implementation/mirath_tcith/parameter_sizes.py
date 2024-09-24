######### Ia
k = 143
m = 16
r = 4
n = 16

rho = 16
mu = 2
q = 16

qbytes = 0.5
qbytes_mu = 1

t = m * r + r * (n - r)
tau = 20
T_OPEN = 113
L = 128 / 8

print(f'\nIa:::(signature size):\t {((t * qbytes + rho * qbytes_mu + 2*L) * tau) + (L * T_OPEN) + (2 * L) + (2 * L) + 8}')
print(f'Ia:::(public key len):\t {(m * n - k) * qbytes + L}')

######### IIIa
k = 288
m = 21
r = 4
n = 21

rho = 24
mu = 2
q = 16

qbytes = 0.5
qbytes_mu = 1

t = m * r + r * (n - r)
tau = 30
T_OPEN = 178
L = 192 / 8

print(f'\nIIIa:::(signature size): {((t * qbytes + rho * qbytes_mu + 2*L) * tau) + (L * T_OPEN) + (2 * L) + (2 * L) + 8}')
print(f'IIIa:::(public key len): {(m * n - k) * qbytes + L}')

######### IIIa
k = 440
m = 25
r = 4
n = 25

rho = 32
mu = 2
q = 16

qbytes = 0.5
qbytes_mu = 1

t = m * r + r * (n - r)
tau = 39
T_OPEN = 247
L = 256 / 8

print(f'\nVa:::(signature size):\t {((t * qbytes + rho * qbytes_mu + 2*L) * tau) + (L * T_OPEN) + (2 * L) + (2 * L) + 8}')
print(f'Va:::(public key len):\t {(m * n - k) * qbytes + L}')