import requests, botConfig

# Replace 'YOUR_BOT_TOKEN' with the token you got from BotFather
BOT_TOKEN = botConfig.TOKEN
# Replace 'YOUR_CHAT_ID' with the chat ID you want to send the message to
CHAT_ID = botConfig.ID

def send_message(text):
    url = f'https://api.telegram.org/bot{BOT_TOKEN}/sendMessage'
    payload = {
        'chat_id': CHAT_ID,
        'text': text
    }
    headers = {
        'Content-Type': 'application/json'
    }
    response = requests.post(url, json=payload, headers=headers)
    return response.json()

# Example usage
message = "Hello, this is a direct message from my code! The run is complete"
response = send_message(message)
print(response)
