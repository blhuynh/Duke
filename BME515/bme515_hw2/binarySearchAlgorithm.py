listOfNumbers = [x*0.1 for x in range(0,10*5)]

def binarySearch(array,value):
	midIndex = math.floor(len(array)/2)
	midValue = array[midIndex]

	if (value == midValue): 
		return midIndex
	elif (value > midValue):
		return binarySearch(array[range(midIndex)],value)
	elif (value < midValue):
		return binarySearch(array[range(midIndex,len(array))],value)
