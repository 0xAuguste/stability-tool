FROM pymesh/pymesh
RUN python -m pip install --upgrade pip
RUN pip3 install matplotlib
RUN pip3 install imageio