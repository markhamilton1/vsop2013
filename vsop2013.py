"""
This module provides basic support for the VSOP2013 ephemerides which can be used to compute the
heliocentric positions and velocities of the 8 planets Mercury, Venus, the Earth-Moon barycenter,
Mars, Jupiter, Saturn, Uranus, Neptune and the dwarf planet Pluto (only in the period from +1500
to +3000).

The ephemerides are provided in 6 sequential text files (ASCII) and can be downloaded at:

    ftp://ftp.imcce.fr/pub/ephem/planets/vsop2013/ephemerides

The following provides the name of each file and the period that it covers:
    VSOP2013.m4000: -4500 to -3000
    VSOP2013.m2000: -3000 to -1500
    VSOP2013.m1000: -1500 to 0
    VSOP2013.p1000: 0 to +1500
    VSOP2013.p2000: +1500 to +3000
    VSOP2013.p4000: +3000 to +4500

While these text files could be used directly, this would have a serious impact on the performance
of computing any significant number of planetary positions. Instead it is suggested that these
files are parsed and stored into binary format for the platform on which they will be used.
NOTE THAT THE BINARY FILE GENERATED FOR ONE PLATFORM MAY NOT BE COMPATIBLE FOR USE ON A DIFFERENT
PLATFORM. Once converted to binary form the text files may be deleted to save storage space where
necessary. The binary form of the files provides a much more efficient mechanism for performing
large numbers of position calculations.

For reference the authors have provided the programs (Fortran) that they used for these purposes.
These are:

VSOP2013_binfile.f :
This program converts the sequential files into direct access files.

VSOP2013_compute.f:
This program computes planetary coordinates from a direct access file.

The file VSOP2013_ctl.txt contains planetary coordinates computed by the program VSOP2013_compute.f
and given as control values for the users.

The file README.pdf also provides this information as well as other info that may be of interest.

This module has been developed to provide a degree of convenience for the Python user, not just
a reference implementation.

TXT_FILES_ROOT is the path to the directory containing the VSOP2013 source text files. You may
change this to reflect where your source text files are actually located. The default is to place
the source text files into a subdirectory called "ephemerides" at the location of the vsop2013.py
script. If the vsop2013.py functionality is being used by other python scripts then this variable
can be changed once the module is loaded.

TXT_FILES is an array of the filenames of the VSOP2013 source text files.

BIN_FILES_ROOT is the path to the directory containing the VSOP2013 binary files. You may
change this to reflect where your binary files are actually located. This variable is similar to
the TXT_FILES_ROOT variable discussed previously.

BIN_FILES is an array of the filenames of the VSOP2013 binary files.

PLANET_NAMES is an array of the names of the supported planets and is used ONLY for the display
of the results by the print_results function.
"""
from array import array
import os
import struct


TXT_FILES_ROOT = "./ephemerides/"
TXT_FILES = ["VSOP2013.m4000",
             "VSOP2013.m2000",
             "VSOP2013.m1000",
             "VSOP2013.p1000",
             "VSOP2013.p2000",
             "VSOP2013.p4000"]

BIN_FILES_ROOT = "./ephemerides/"
BIN_FILES = ["VSOP2013.m4000.pybin",
             "VSOP2013.m2000.pybin",
             "VSOP2013.m1000.pybin",
             "VSOP2013.p1000.pybin",
             "VSOP2013.p2000.pybin",
             "VSOP2013.p4000.pybin"]


PLANET_NAMES = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

def print_results(iplanet, jde, r):
    name = PLANET_NAMES[iplanet]
    print("{:8s} JD{:9.1f}  X : {: 17.13f} ua    Y : {: 17.13f} ua    Z : {: 17.13f} ua".format(name, jde, r[0], r[1], r[2]))
    print("                      X': {: 17.13f} ua/d  Y': {: 17.13f} ua/d  Z': {: 17.13f} ua/d".format(r[3], r[4], r[5]))


class VSOP2013File:

    N_PLANETS = 9
    DATE_RANGES = [(77294.5, 625198.5),    # -4500 to -3000
                  (625198.5, 1173102.5),   # -3000 to -1500
                  (1173102.5, 1721006.5),  # -1500 to 0
                  (1721006.5, 2268910.5),  # 0 to +1500
                  (2268910.5, 2816814.5),  # +1500 to +3000
                  (2816814.5, 3364718.5)]  # +3000 to +4500

    def __init__(self):
        self.bfref = None
        self.bfi = None
        self.doff = None
        self.idf = None
        self.t1 = None
        self.t2 = None
        self.delta = None
        self.nintv = None
        self.ncoef = None
        self.loc = None

    # text file support methods
    #
    # These methods are here to support the conversion of the coefficient text files to binary files.
    # This conversions decreases the size of the coefficient files, eliminates the text parsing by putting
    # the data into the proper format, and makes random access within a file possible, dramatically
    # improving performance.

    def txtfile_path(self, i):
        """
        Build a source text file path.
        :param i: index of the source text file in the TXT_FILES array
        :return: the source text file path (None if error)
        """
        if (i >= 0) and (i < self.N_PLANETS):
            return os.path.join(TXT_FILES_ROOT, TXT_FILES[i])
        return None

    def txtfile_exists(self, i):
        """
        Test if a source text file exists.
        :param i: index of the source text file in the TXT_FILES array
        :return: true=exists and is a file
        """
        tfp = self.txtfile_path(i)
        if tfp is not None:
            return os.path.isfile(tfp)
        return False

    def bin_txtfile(self, i):
        """
        Convert a source text file to a binary file.
        :param i: index of the source text file in the TXT_FILES array
        """
        self.bfref = None
        self.bfi = None
        self.doff = None
        self.idf = None
        self.t1 = None
        self.t2 = None
        self.delta = None
        self.nintv = None
        self.ncoef = None
        self.loc = None
        tfp = self.txtfile_path(i)
        if tfp is None:
            raise ValueError("Unable to build text file path with i = {}!".format(i))
        bfp = self.binfile_path(i)
        if bfp is None:
            raise ValueError("Unable to build binary file path with i = {}!".format(i))
        if not os.path.isfile(tfp):
            raise ValueError("File '{}' not found!".format(tfp))
        print("Creating '{}' from '{}'".format(bfp, tfp))
        with open(tfp, 'r') as tfref:
            self.__load_idf(tfref)
            self.__load_t1(tfref)
            self.__load_t2(tfref)
            self.__load_delta(tfref)
            self.__load_nintv(tfref)
            self.__load_ncoef(tfref)
            self.__load_loc(tfref)
            self.bfref = open(bfp, 'wb')
            self.__write_binfile_header()
            self.__bin_records(tfref)
            self.bfref.close()
            self.bfref = None

    def __load_idf(self, tfref):
        line = tfref.readline()
        try:
            v = int(line.strip())
        except:
            raise ValueError("Invalid value for idf!")
        if v != 2013:
            raise ValueError("Invalid value for idf!")
        self.idf = "vsop2013"

    def __load_t1(self, tfref):
        line = tfref.readline()
        try:
            self.t1 = float(line.strip())
        except:
            raise ValueError("Invalid value for t1!")

    def __load_t2(self, tfref):
        line = tfref.readline()
        try:
            self.t2 = float(line.strip())
        except:
            raise ValueError("Invalid value for t2!")

    def __load_delta(self, tfref):
        line = tfref.readline()
        try:
            self.delta = float(line.strip())
        except:
            raise ValueError("Invalid value for delta!")

    def __load_nintv(self, tfref):
        line = tfref.readline()
        try:
            self.nintv = int(line.strip())
        except:
            raise ValueError("Invalid value for nintv!")

    def __load_ncoef(self, tfref):
        line = tfref.readline()
        try:
            self.ncoef = int(line.strip())
        except:
            raise ValueError("Invalid value for ncoef!")

    def __load_loc(self, tfref):
        line = tfref.readline()
        p1 = line.strip().split(" ")
        line = tfref.readline()
        p2 = line.strip().split(" ")
        line = tfref.readline()
        p3 = line.strip().split(" ")
        self.loc = [[], [], []]
        for i in range(0, len(p1)):
            self.loc[0].append(int(p1[i]))
            self.loc[1].append(int(p2[i]))
            self.loc[2].append(int(p3[i]))

    def __write_binfile_header(self):
        self.bfref.write(struct.pack('8s', self.idf))
        self.bfref.write(struct.pack('fffii', self.t1, self.t2, self.delta, self.nintv, self.ncoef))
        self.bfref.write(struct.pack('9i', *self.loc[0]))
        self.bfref.write(struct.pack('9i', *self.loc[1]))
        self.bfref.write(struct.pack('9i', *self.loc[2]))

    def __bin_record(self, tfref):
        line = tfref.readline()
        p = line.strip().split(" ")
        try:
            d1 = float(p[0])
            d2 = float(p[-1])
        except:
            raise ValueError("Invalid value for date!")
        n = self.ncoef
        coefs = []
        while n > 0:
            line = tfref.readline()
            line = line.strip().replace("  ", " ")
            p = line.split(" ")
            i = 0
            while i < len(p):
                try:
                    v = float(p[i])
                    e = float(p[i+1])
                    c = v * (10 ** e)
                    coefs.append(c)
                    n -= 1
                    i += 2
                except:
                    raise ValueError("Invalid coefficient!")
        self.bfref.write(struct.pack('2f', d1, d2))
        self.bfref.write(struct.pack('%id' % len(coefs), *coefs))

    def __bin_records(self, tfref):
        for i in range(0, self.nintv):
            if i % 1000 == 0:
                print("{} / {}".format(i, self.nintv))
            self.__bin_record(tfref)

    # bin file support methods

    def binfile_path(self, i):
        """
        Build a binary file path.
        :param i: index of the file name in the BIN_FILES array
        :return: the bin file path (None if error)
        """
        if (i >= 0) and (i < self.N_PLANETS):
            return os.path.join(BIN_FILES_ROOT, BIN_FILES[i])
        return None

    def binfile_exists(self, i):
        """
        Test if a binary file exists.
        :param i: index of the binary file name in the BIN_FILES array
        :return: true if the file exists
        """
        bfp = self.binfile_path(i)
        if bfp is not None:
            return os.path.isfile(bfp)
        return False

    def close_binfile(self):
        """
        Close the current binary file if necessary.
        """
        if self.bfref is not None:
            self.bfref.close()
            self.bfref = None
            self.bfi = None

    def open_binfile(self, i):
        """
        Open the specified binary file.
        :param i: index of the binary file in the BIN_FILES array
        """
        self.close_binfile()
        bfp = self.binfile_path(i)
        self.bfref = open(bfp, 'rb')
        self.bfi = i
        self.__read_binfile_header()
        self.doff = self.bfref.tell()

    def __read_binfile_header(self):
        self.idf = struct.unpack('8s', self.bfref.read(8))[0]
        if self.idf != "vsop2013":
            raise ValueError("Invalid file id! {}".format(self.idf))
        hdr = struct.unpack('fffii', self.bfref.read(5 * 4))
        self.t1 = hdr[0]
        self.t2 = hdr[1]
        self.delta = hdr[2]
        self.nintv = hdr[3]
        self.ncoef = hdr[4]
        self.loc = [None, None, None]
        self.loc[0] = struct.unpack('9i', self.bfref.read(9 * 4))
        self.loc[1] = struct.unpack('9i', self.bfref.read(9 * 4))
        self.loc[2] = struct.unpack('9i', self.bfref.read(9 * 4))

    def calculate_for(self, ip, jde):
        """
        Attempt to calculate the position and velocity for the specified planet.
        :param ip: index of planet to perform calculation for
        :param jde: the julian date to perform the calculation for
        :return: an array of the X, Y, and Z heliocentric position of the planet and the X', Y', and Z' heliocentric velocity
        """
        if (ip < 0) or (ip >= self.N_PLANETS):
            raise ValueError("Invalid planet index! Is outside range of [{}:{}]".format(0, self.N_PLANETS))
        i = None
        rng = None
        for i, r in enumerate(self.DATE_RANGES):
            if (jde >= r[0]) and (jde < r[1]):
                rng = r
                break
        if rng is None:
            raise ValueError("Invalid jde! Is outside range of [{}:{}]".format(self.DATE_RANGES[0][0], self.DATE_RANGES[-1][1]))
        # open bin file and read header
        if i != self.bfi:
            self.open_binfile(i)
        # calculate offsets
        irec = int((jde - self.t1) / self.delta)
        coefs_size = self.ncoef * 8
        off = irec * ((2 * 4) + coefs_size)
        # read the record for the date
        self.bfref.seek(self.doff + off)
        dj1, dj2 = struct.unpack('2f', self.bfref.read(2 * 4))
        coefs = array('d')
        coefs.fromfile(self.bfref, coefs_size)
        iad = self.loc[0][ip] - 1
        ncf = self.loc[1][ip]
        nsi = self.loc[2][ip]
        delta2 = self.delta / nsi
        ik = int((jde - dj1) / delta2)
        if ik == nsi:
            ik -= 1
        iloc = iad + (6 * ncf * ik)
        dj0 = dj1 + (ik * delta2)
        x = (2.0 * (jde - dj0) / delta2) - 1.0
        tn = [0.0] * ncf
        tn[0] = 1.0
        tn[1] = x
        for i in range(2, ncf):
            tn[i] = (2.0 * x * tn[i - 1]) - tn[i - 2]
        r = [0.0] * 6
        for i in range(0, 6):
            for j in range(0, ncf):
                jp = ncf - j - 1
                jt = iloc + (ncf * i) + jp
                r[i] = r[i] + (tn[jp] * coefs[jt])
        return r


if __name__ == "__main__":


    vsop2013 = VSOP2013File()

    # generate bin files from source text files if necessary
    for i in range(0, len(TXT_FILES)):
        if vsop2013.txtfile_exists(i) and not vsop2013.binfile_exists(i):
            vsop2013.bin_txtfile(i)

    # calculate positions and velocities as a control set for algorithm verification
    ndat = 5
    YEAR = (-4500, -3000, -1500, 0, 1500, 3000)
    TZERO = (77432.5, 625307.5, 1173182.5, 1721057.5, 2268932.5, 2816818.5)
    step = 136798.0
    for i, tzero in enumerate(TZERO):
        print
        print("*** {:>5d} to {:<5d} ***".format(YEAR[i], YEAR[i] + 1500))
        print
        for ip in range(0, vsop2013.N_PLANETS):
            for n in range(0, ndat):
                jd = tzero + (n * step)
                r = vsop2013.calculate_for(ip, jd)
                if r[0] != 0.0:
                    print_results(ip, jd, r)
    vsop2013.close_binfile()

