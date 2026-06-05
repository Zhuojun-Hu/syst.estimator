Step 1: Start a jupyter lab server at port 9527 (or what ever you like)
```
source start_jupyter_server.csh
```
You should see an output like
```
http://localhost:9527/lab?token=XXXXXXXXXXX
```
---
Step 2: Make a new connection to sukap, check out port 9527. For example, if you use VS code, make a new connection:
```
Ctrl+Alt+P
```
then
```
Remote-SSH: Connect to Host...
ssh username@sukap02
```
---
Step 3: From ports, find 
```
localost:9527
```
from Forwarded Address.
Then select 
```
Open in browser
```
Enter the token from step 1. You should be able to open a jupyter lab on VS code browser.