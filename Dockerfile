FROM python:3.12-slim
RUN useradd --create-home --shell /bin/bash app_user
WORKDIR /home/app_user
COPY . .
RUN pip install --no-cache-dir .
USER app_user
CMD ["bash"]