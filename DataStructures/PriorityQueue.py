# code taken from geeksforgeeks
# https://www.geeksforgeeks.org/priority-queue-in-python/

# A simple implementation of Priority Queue
# using Queue.
class PriorityQueue(object):
    def __init__(self):
        self.queue = []
  
    def __str__(self):
        return ' '.join([str(i) for i in self.queue])
  
    # for checking if the queue is empty
    def isEmpty(self):
        return len(self.queue) == 0
  
    # for inserting an element in the queue
    def insert(self, data):
        self.queue.append(data)
  
    # for popping an element based on Priority
    def pop(self):
        try:
            min = 0
            for i in range(len(self.queue)):
                if self.queue[i].totalCost < self.queue[min].totalCost:
                    min = i
            item = self.queue[min]
            del self.queue[min]
            return item
        except IndexError:
            print()
            exit()


if __name__ == '__main__':
    myQueue = PriorityQueue()
    myQueue.insert(12)
    myQueue.insert(1)
    myQueue.insert(14)
    myQueue.insert(7)
    print(myQueue)            
    while not myQueue.isEmpty():
        print(myQueue.pop()) 