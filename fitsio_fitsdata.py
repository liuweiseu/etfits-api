import sys
import FitsData
import fitsio

def pure_python_test(filename):
  data = FitsData.FitsData(filename)

if __name__ == '__main__':
  filename = sys.argv[1]
  pure_python_test(filename)
