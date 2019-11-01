from orbittools import *

def main():
    print('Expected value for julDay(2000,1,1,12,0,0) is 2.451545e6')
    print('%.15e'%julDay(2000,1,1,12,0,0))
    print('Expected value for julDay(2013,6,4,15,9,23) is 2.456448131516204e6')
    print('%.15e'%julDay(2013,6,4,15,9,23))
    print('Expected value for modJulDay(1996,12,31,0,0,0) is 50448.0')
    print(modJulDay(1996,12,31,0,0,0))
    print('Expected value for toModJulDay(toJulDay(modJulDay(2015,7,14,12,0,0))) is 57217.5')
    print(toModJulDay(toJulDay(modJulDay(2015,7,14,12,0,0))))
    
    print('Additional testing not in notebook')
    print('Expected value for toJulDay(modJulDay(2015,7,14,12,0,0))) is 2457218.0')
    print(toJulDay(modJulDay(2015,7,14,12,0,0)))
    print('Test of centT: Centuries past J2000')
    print(centT(toJulDay(modJulDay(2315,7,14,12,0,0))))


if __name__ == "__main__":
    main()
