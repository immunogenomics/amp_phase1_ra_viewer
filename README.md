ampviewer
=========

Check it out here: https://amp.quickgene.net

Contact [Kamil] for the username and password.

[Kamil]: mailto:kslowikowski@gmail.com?subject=ampviewer

Overview
--------

ampviewer is an app for viewing AMP data.

<img height="300px" src="screenshot.png" />

Deploy to our web server
------------------------

Run the app on your laptop when you're developing new features.

When you're ready to launch the app on the public website, copy it to our Amazon EC2 server:

```bash
src="/PATH/TO/ampviewer/"
des="ig:/srv/shiny-server/ampviewer/"

rsync -avh \
    --exclude=*.sh \
    --exclude=*.Rproj \
    --exclude=data/*.rds \
    --exclude=.* \
    "$src" "$des"
```

Get ssh access to our web server
--------------------------------

1. Ask Kamil to add your username to the server.
2. Create a private/public key pair on your laptop.
3. Share the public key with Kamil.

Create a new private/public key pair like this:

```
ssh-keygen -t rsa
```

You will see the output shown below. Please follow these instructions:

1. Please change `/Users/YOUR_USERNAME/.ssh/id_rsa` to `/Users/YOUR_USERNAME/.ssh/immunogenomics_id_rsa`

2. You can use an empty passphrase if you wish (just press <kbd>Enter</kbd>).

3. This creates 2 files:

    - a **private** key (never share with anyone) `~/.ssh/immunogenomics_id_rsa`
    
    - a **public** key (you can share freely with the world) `~/.ssh/immunogenomics_id_rsa.pub`

4. Please share the public key with Kamil so he can give you access to the immunogenomics.io server.

```
Generating public/private rsa key pair.
Enter file in which to save the key (/Users/YOUR_USERNAME/.ssh/id_rsa): /Users/YOUR_USERNAME/.ssh/immunogenomics_id_rsa
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /Users/YOUR_USERNAME/.ssh/immunogenomics_id_rsa.
Your public key has been saved in /Users/YOUR_USERNAME/.ssh/immunogenomics_id_rsa.pub.
The key fingerprint is:
SHA256:d+DRaMtQEyECvQBna7R7BRihuvfE2Y7wos2u+/GX2LM YOUR_USERNAME@hostname
The key's randomart image is:
+---[RSA 2048]----+
|  ..B*o . =o     |
|   =oo.o o +     |
|  . +. .o = .    |
| . . ... = +     |
|.   . . S = .    |
| . . +   . .     |
|. + +o..         |
| +.B.o=          |
|+*=.=.Eo         |
+----[SHA256]-----+
```

On your laptop, please add this to your `~/.ssh/config` file:

```
Host ig
    PubkeyAuthentication yes
    User USERNAME
    IdentityFile ~/.ssh/immunogenomics_id_rsa
    HostName 34.196.99.74
```

On the immunogenomics.io server, let's add you as a new user and give you ssh access:

```bash
# Add a user
sudo adduser USERNAME

# Change to that user
sudo su USERNAME

# Setup ssh
cd
mkdir .ssh
touch .ssh/authorized_keys
sudo chmod 0700 .ssh
sudo chmod 0600 .ssh/authorized_keys
echo 'ssh-rsa AAAAB3NzaC1yc2EA...' >> ~/.ssh/authorized_keys
```

Finally, let's add the new user to the `shiny` group:

```bash
sudo usermod -a -G shiny USERNAME
```


