import os
import requests
import logging
from requests.exceptions import RequestException
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO)

# Load sensitive data from environment variables
BOT_TOKEN = os.getenv('TOKEN')
CHAT_ID = os.getenv('ID')

def send_message(text, parse_mode=None, disable_notification=False):
    url = f'https://api.telegram.org/bot{BOT_TOKEN}/sendMessage'
    payload = {
        'chat_id': CHAT_ID,
        'text': text
    }
    headers = {
        'Content-Type': 'application/json'
    }
    
    try:
        response = requests.post(url, json=payload, headers=headers)
        response.raise_for_status()
        logging.info(f"Message sent: {text}")
        return response.json()
    except RequestException as e:
        logging.error(f"Failed to send message: {e}")
        return None

# Example usage
if __name__ == "__main__":
    message = "Hello, this is a direct message from my code! The run is complete"
    response = send_message(message)
    print(response)
