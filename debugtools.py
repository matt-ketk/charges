import time

'''
just a module used to compare runtimes of functions. will probably be used in later episodes of optimization. -matt 
'''

# millis = int(round(time.time() * 1000))

def getFuncRuntime(func, *args):
  t = getCurrentTimeMillis()
  func(*args)
  return getCurrentTimeMillis() - t

def getCurrentTimeMillis():
  return int(round(time.time() * 1000))