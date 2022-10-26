import unittest
import colonyMaskCreation as s
import numpy as np
from numpy import testing

class ColonyMask(unittest.TestCase):
    def testSize_cell_matrix_record(self):
        self.assertEqual(np.shape(s.cell_matrix_record),(s.J,s.J,s.N), 'Expected to be JJN' )

    def test_check_neighbours(self):
        cell_matrix = np.ones((5,5))
        y = 2
        x = 3
        result = s.check_neighbours(cell_matrix, x, y)
        
        expected = np.array([[1,1,1],[1,np.nan,1],[1,1,1]])
        testing.assert_array_equal(result, expected)
if __name__ =='__main__':
    unittest.main()

    