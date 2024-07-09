# ex-fem-user

Explicit Finite Element Method for Simulating Dynamic Problems

## Getting Started
To Run as a User:
- Create a virtual environment
- Clone the repo
- Install the requirements
- Run!
<pre>
mkdir test-ex-fem
cd test-ex-fem
python3.8 -m venv env
env\Scripts\activate.bat
git clone https://github.com/kinfungchan/ex-fem-user.git
pip install -r ex-fem-user\requirements.txt
python ex-fem-user\ex-fem\main.py
</pre>
Now you can run a simple benchmark where a wave propagates through a 2D Mesh