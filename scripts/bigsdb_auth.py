#!/usr/bin/env python3
"""
bigsdb_auth.py - one-time helper to get a BIGSdb-Pasteur OAuth access token.

Run this ONCE on your own computer (not in CI):

    pip install requests requests-oauthlib
    python scripts/bigsdb_auth.py

You need the consumer key and secret that the BIGSdb-Pasteur team emailed
you. The script prints a link; open it, log in to BIGSdb-Pasteur, click
Authorize, and paste the verification code back here. It then prints the
access token + secret to store as GitHub secrets. The access token does not
expire (unless you revoke it in your BIGSdb account); the updater uses it to
obtain short-lived session tokens on each run.

Nothing is written to disk.
"""
import argparse
import getpass
import sys

DB = "pubmlst_diphtheria_seqdef"
API = f"https://bigsdb.pasteur.fr/api/db/{DB}"
WEB = "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--api", default=API)
    p.add_argument("--web", default=WEB, help="BIGSdb web script URL")
    p.add_argument("--db", default=DB)
    args = p.parse_args()

    try:
        from requests_oauthlib import OAuth1Session
    except ImportError:
        sys.exit("Install first: pip install requests requests-oauthlib")

    key = input("Consumer key (client id): ").strip()
    secret = getpass.getpass("Consumer secret (client secret, hidden): ").strip()

    # BIGSdb reads the OAuth fields from the query string, not the
    # Authorization header ("No consumer key submitted" otherwise).
    # 1. request token
    s = OAuth1Session(key, client_secret=secret, callback_uri="oob",
                      signature_type="QUERY")
    r = s.get(f"{args.api}/oauth/get_request_token", timeout=60)
    if r.status_code != 200:
        sys.exit(f"Request token failed ({r.status_code}): {r.text[:300]}\n"
                 "Check the consumer key/secret.")
    req = r.json()

    # 2. user authorises in browser
    url = (f"{args.web}?db={args.db}&page=authorizeClient"
           f"&oauth_token={req['oauth_token']}")
    print("\nOpen this link, log in, and click Authorize:\n\n  " + url + "\n")
    verifier = input("Verification code shown after authorising: ").strip()

    # 3. access token
    s = OAuth1Session(key, client_secret=secret,
                      resource_owner_key=req["oauth_token"],
                      resource_owner_secret=req["oauth_token_secret"],
                      verifier=verifier, signature_type="QUERY")
    r = s.get(f"{args.api}/oauth/get_access_token", timeout=60)
    if r.status_code != 200:
        sys.exit(f"Access token failed ({r.status_code}): {r.text[:300]}\n"
                 "The code may have expired - run again and paste it quickly.")
    acc = r.json()

    # 4. sanity check: can we get a session token?
    s = OAuth1Session(key, client_secret=secret,
                      resource_owner_key=acc["oauth_token"],
                      resource_owner_secret=acc["oauth_token_secret"],
                      signature_type="QUERY")
    ok = s.get(f"{args.api}/oauth/get_session_token", timeout=60).status_code == 200

    print("\nSuccess." if ok else "\nGot an access token, but the session-token "
          "check failed - check you are registered to the database.")
    print("Add these as GitHub repository secrets "
          "(Settings > Secrets and variables > Actions):\n")
    print(f"  BIGSDB_CONSUMER_KEY     = {key}")
    print("  BIGSDB_CONSUMER_SECRET  = (the consumer secret you entered)")
    print(f"  BIGSDB_ACCESS_TOKEN     = {acc['oauth_token']}")
    print(f"  BIGSDB_ACCESS_SECRET    = {acc['oauth_token_secret']}")
    print("\nTreat these like passwords. Clear your terminal afterwards.")


if __name__ == "__main__":
    main()
