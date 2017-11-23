import h5py


def text_zwrt_quanty_inp(num_b, num_e, h1e, v2e):
    with open("quanty.inp", "w") as f:
        f.write("{}\n".format(num_b))
        f.write("{}\n".format(num_e))
        for h1 in h1e:
            line = ""
            for h12 in h1:
                line += "{:18.12} {:18.12}".format(h12.real, h12.imag)
            f.write(line+"\n")
        for v1 in v2e:
            for v12 in v1:
                for v123 in v12:
                    line = ""
                    for v1234 in v123:
                        line += "{:18.12} {:18.12}".format( \
                                v1234.real, v1234.imag)
                    f.write(line+"\n")


def get_h1e_v2e(fname="EMBED_HAMIL_1.h5"):
    with h5py.File(fname, "r") as f:
        h1e = f["/H1E"][()].T
        v2e = f["/V2E"][()].T
    return h1e, v2e


if __name__ == "__main__":
    h1e, v2e = get_h1e_v2e()
    num_b = h1e.shape[0]
    num_e = 6
    text_zwrt_quanty_inp(num_b, num_e, h1e, v2e)
