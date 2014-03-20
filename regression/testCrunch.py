import os
import sys
import string
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import subprocess
import argparse
import datetime

class TecplotFile:
    """
    Tecplot file class
    """
    def __init__(self,name,path):
        self.name = name
        self.path = path
        self.pathname = os.path.join(path,name)
        self.File = None
        self.lines = []
        
    def __str__(self):
        words = 'Tecplot file: {0}'.format(self.pathname)
        return words

    def _readlines(self):

        self.File = open(self.pathname,'r')
        self.lines = self.File.readlines()
        self.File.close()
    
    def get_values(self,target_vars):

      self._readlines() 
        
      values = [ [] for var in target_vars ]
      target = None

      for line in self.lines:

          if line.lstrip()[0] =='#':        
        
              if 'Units:' in line:
                  units = line.split(':')[1].lstrip()
              elif 'Time series at grid cell:': 
                  cell  = line.split(':')[1].lstrip()
                  cellx,celly,cellz = cell.split()
    
          elif 'TITLE' in line:
          
              title = line.strip().split('=')[1].translate(None,'"').lstrip()
    
          elif 'VARIABLES' in line:
        
              allvariables = line.strip().split('=')[1].split(',')
              allvariables = [var.translate(None,'"').lstrip().strip() for var in allvariables]
              allvariables = [var for var in allvariables if var != '']

          elif 'ZONE F=POINT' in line:
        
              ipoints = int(line.strip().split(',')[1].split('=')[1])
              if len(line.strip().split(',')) > 2:
                  jpoints = int(line.strip().split(',')[2].split('=')[1])
                  if len(line.strip().split(',')) > 3:
                      kpoints = int(line.strip().split(',')[3].split('=')[1])

          else:
        
              if len(allvariables) == 0:
                  print "no variables provided in file"
              elif target == None:
                  target = [allvariables.index(var) for var in target_vars]

              for i,j in enumerate(target):          
                  values[i] = values[i] + [float(line.split()[j])]
  
      return values

# end TecplotFile class

class CrunchTest:
    """
    Test class
    """    
    def __init__(self,name,path,executable,pathexe):
        self.name = name
        self.path = path 
        self.executable = executable
        self.pathexe = pathexe
        self.subtests = []

    def __str__(self):
        words =  ' Test name: {0}'.format(self.name)
        words += '\n Test location: {0}'.format(self.path)
        words += '\n Executable: {0}'.format(self.executable)
        words += '\n Executable location: {0}'.format(self.pathexe)
        
        return words

    def _read_params(self):

        paramFileName = self.name + '.par'
        paramFileName = os.path.join(self.path,paramFileName)

        File = open(paramFileName,'r')
        self.lines = File.readlines()
        File.close()

        for line in self.lines:

          if 'testcf' in line:

              filename   = line.split('=')[1].lstrip().rstrip()

              if ',' in line.split('=')[2]:               
                print line.split('=')[2]
                print 'here' 
                variables  = line.split('=')[2].split(',')
                tolerances = line.split('=')[3].split(',')
              else:
                print line.split('=')[2]
                print 'there'
                variables  = [line.split('=')[2]]
                tolerances = [line.split('=')[3]]
              benchmark  = line.split('=')[4].lstrip().rstrip()
              typetest  =  line.split('=')[5].split(',') 
                    
              # if len(variables) > len(tolerances):
              
              variables = [var.lstrip().rstrip() for var in variables]
              tolerances = [float(tol) for tol in tolerances]
              typetest = [tt.lstrip().rstrip() for tt in typetest]
                
              self.subtests = self.subtests + [[filename,variables,tolerances,benchmark,typetest]]
        
        print self.subtests
                
    def _run_CrunchFlow(self):
                       
        CrunchFlow_exe = [os.path.join(self.pathexe,self.executable)]

        err = subprocess.call(CrunchFlow_exe)
        
        return err
        
    def _compare(self):

        print 'read params'

        self._read_params()

        failed = 0

        for subtest in self.subtests:

           print subtest
                  
           filename  = subtest[0]
           variables = subtest[1]
           tolerances = subtest[2]
           benchmark  = subtest[3] 
           typetest   = subtest[4]
        
           calc = TecplotFile(filename,self.path)           
           print 'tecplot file output',filename
           calc = np.array(calc.get_values(variables))
          
           print 'past calc read'

           gold = TecplotFile(benchmark,self.path)
           print 'tecplot file benchmark',filename
           gold = np.array(gold.get_values(variables))

           print 'past gold read'

           if typetest[1] == 'time':

                timevar = np.split(calc,len(variables),axis=0)[0][0]
                timegold = np.split(gold,len(variables),axis=0)[0][0]
                
                newcalc = [timegold]
                
                for i in range(len(variables))[1:]:

                   yval = np.split(calc,len(variables),axis=0)[i][0]

                   s = InterpolatedUnivariateSpline(timevar, yval)
                   newcalc = np.append(newcalc,[s(timegold)],axis=0)
                
                calc = newcalc

           print len(calc), len(gold)
           diff = calc - gold

           print typetest[0]

           if typetest[0] == 'relative':
                diff = diff / gold
                
           #if  np.linalg.norm(diff,axis=1) < np.array(subtest[2]).all():
           print 'past interpolation'

           for i,var in enumerate(variables):
                    
              norm = np.linalg.norm(np.split(diff,len(variables),axis=0)[i])

              if  norm < tolerances[i]:
                    print 'test passed for {0}'.format(var)              
                    print 'norm {0}\n'.format(norm)  

              else:
                    print 'test failed for {0}'.format(var)
                    print 'norm {0}\n'.format(norm)
                    failed = failed + 1

           print 'past norm'
        
        return failed

    def _clean(self):
          
        files = os.listdir(os.getcwd())

        for filename in files:
          if ".out" in filename:
             os.remove(filename)

        return

    def run(self):
                     
        os.chdir(self.path)
        print "running example in {0}".format(os.getcwd())
               
        err = self._run_CrunchFlow()

        failed = self._compare()

        self._clean()

        os.chdir('..')
    
        return failed

# end CrunchTest class

class CrunchTestSuite:
    """
    CrunchTestSuite class
    """    
    def __init__(self,name):
         self.name = name
         self.executable = ''
         self.pathexe = ''
         self.tests = []
         self.failed = []

   ## def _read_params(self):
    
         paramFileName = self.name + '.par'
	 paramFileName = os.path.join(paramFileName)

         File = open(paramFileName,'r')
	 lines = File.readlines()
	 File.close()

	 for line in lines:
	
                  if 'executablecf' in line:
                  
                      self.executable  =  line.split('=')[1].split(',')[0].lstrip().rstrip()
                      
                      pathexe     =  line.split('=')[1].split(',')[1].lstrip().rstrip()
                      if pathexe == 'root':
                         self.pathexe = '' 
                      else:
                         self.pathexe = pathexe

	          elif 'testcf' in line:
	
	              name  =  line.split('=')[1].split(',')[0].lstrip().rstrip()
                      path  =  line.split('=')[1].split(',')[1].lstrip().rstrip()

	              self.tests = self.tests + [[name,path]]

    def run(self):

         current = os.getcwd()


         # go to root to chdir into path to executable	 
	 if 'regression' in current:
	     os.chdir('..')
	     current = os.getcwd()

         pathexe = os.path.join(current,self.pathexe)

         # go back to regression dir, assume examples there        
         os.chdir('regression')
         current = os.getcwd()

         for test in self.tests:

             name    = test[0]
             path    = os.path.join(current,test[1])

             failed = CrunchTest(name,path,self.executable,pathexe).run()
             self.failed = self.failed + [failed] 
         
    def report(self):

         now = datetime.datetime.now()
         now = now.strftime("%Y%m%d_%H%M%S")

         repfile = open(self.name+'_'+now+'.report','w')
         bar = '-'*40+'\n'

         repfile.write(bar)
         repfile.write('        Test suite summary \n')
         repfile.write(bar)

         for (test,status) in zip(self.tests,self.failed):

             minireport = 'test: {0}, fails: {1} \n'.format(test[0],status)
             repfile.write(minireport)

         repfile.write(bar)
         repfile.close()

if __name__ == "__main__":

        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='''
command line launch of crunch regression tests.''',
                                     epilog='''
only name of test suite supported for now

    - run test suite for Crunch "NAME":
        ./testCrunch.py --test_suite NAME

''')    
        parser.add_argument('-t', '--test_suite', metavar='testsuitename', nargs=1,
                            required=False, default=None,
                            help='name of test suite, no file extension')

        options = parser.parse_args()
      
        if len(options.test_suite):
          testsuitename = options.test_suite[0]
        else:
          testsuitename = 'all'

        try:
            
            a = CrunchTestSuite(testsuitename)
            a.run()
            a.report()

        except:

            print "this whole thing did not go as planned"
