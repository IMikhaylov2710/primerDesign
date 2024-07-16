import threading

def checkCredentials(task, quota):
    threading.Timer(1.0, checkCredentials).start()
    if task >= quota:
        print('No more nodes to allocate')
    else:
        print(f'Allocated {task} out of {quota} nodes')

checkCredentials()
