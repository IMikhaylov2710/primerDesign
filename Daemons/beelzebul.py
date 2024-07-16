import threading

def checkCredentials(task, quota):
    threading.Timer(1.0, checkCredentials).start()
    task += 1
    if task >= quota:
        print('No more nodes to allocate')
    else:
        print(f'Allocated {task} out of {quota} nodes')

checkCredentials(0, 10)
