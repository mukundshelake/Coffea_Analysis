import time

def long_running_task():
    print("Task started...")
    time.sleep(10)  # Sleep for 25 seconds
    print("Task completed!")

if __name__ == "__main__":
    long_running_task()