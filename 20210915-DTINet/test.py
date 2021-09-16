import socket
client = socket.socket()
client.connect(('172.17.0.2', 8777))
client.send(b'hello')
client.close()